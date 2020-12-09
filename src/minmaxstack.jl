struct MinMaxSpec{T,IM}
    times::T
    interpmode::IM
end

interpmode(mm::MinMaxSpec) = mm.interpmode

function chooseinterparrays(spec::MinMaxSpec, step, tfrac) 
    # Get layers, indices and interpolation fraction fo this specific
    # fraction of the timestep
    minfrac, maxfrac = map(t -> Second(t) / Second(Day(1)), spec.times)
    # The last section, larger than the max time
    if tfrac >= maxfrac
        stackkeys = reverse(keys(spec.times))
        window = maxfrac, minfrac + 1
        inds = 1, 2 
    elseif tfrac < minfrac
        stackkeys = reverse(keys(spec.times))
        window = maxfrac - 1, minfrac
        inds = 0, 1 
    else
        stackkeys = keys(spec.times)
        window = minfrac, maxfrac
        inds = 1, 1 
    end
    wtfrac = -(window[1] - tfrac) / (window[2] - window[1]) 
    NamedTuple{stackkeys}(inds), wtfrac
end

struct MinMaxFracs{I,F,IM}
    indices::I
    frac::F
    interpmode::IM
end
MinMaxFracs(spec, step, trfrac) = begin
    inds, frac = chooseinterparrays(spec, step, trfrac)
    imode = spec.interpmode
    MinMaxFracs{typeof(inds),typeof(frac),typeof(imode)}(inds, frac, imode)
end

frac(mm::MinMaxFracs) = mm.frac
indices(mm::MinMaxFracs) = mm.indices
interpmode(mm::MinMaxFracs) = mm.interpmode

"""
    MimMaxInterpStack

Provides interpolation for a value (usually temperature) that usually
has minimum and maximum values provided in separate layers of the
stack. Other layers in the stack are interpolated, or not,
as for `InterpStack`.
"""
struct MinMaxInterpStack{S,I,MMS<:NamedTuple} <: AbstractGeoStack{S}
    stacks::S
    interpmodes::I
    frac::Float64
    specs::MMS
end
MinMaxInterpStack(stacks, interpmodes, tfrac, specs::NamedTuple{<:Any,<:Tuple{Vararg{<:MinMaxSpec}}}, timestep) = begin
    minmaxfracs = map(spec -> MinMaxFracs(spec, timestep, tfrac), specs)
    interpseries[i] = MinMaxInterpStack(stacks, interpmodes, frac, minmaxfracs)
end

stacks(s::MinMaxInterpStack) = s.stacks
interpmodes(s::MinMaxInterpStack) = s.interpmodes
frac(s::MinMaxInterpStack) = s.frac
specs(s::MinMaxInterpStack) = s.specs

GeoData.refdims(s::MinMaxInterpStack) = refdims(first(stacks(s)))
GeoData.metadata(s::MinMaxInterpStack) = metadata(first(stacks(s)))
GeoData.window(s::MinMaxInterpStack) = GeoData.window(first(stacks(s)))
GeoData.childtype(s::MinMaxInterpStack) = GeoData.childtype(first(stacks(s)))


Base.getindex(s::MinMaxInterpStack, key::Symbol) =
    if key in keys(s.specs)
        spec = s.specs[key]
        frac_ = frac(spec)
        intmode = interpmode(spec)
        arrays = map((k, i) -> s.stacks[i][k], keys(spec.indices), spec.indices)
        imode = interpmode(spec)
        if frac_ < 1
            frac_ += 1
        end
        valweights = Interpolations.value_weights(imode.degree, frac_)
        weights = Interpolations.WeightedIndex(axes(arrays, 1), valweights)
        InterpArray(arrays, imode, frac_, interpolate(arrays, imode), weights)
    else
        interparray(s, key)
    end

# TODO include non-minmax keys
Base.keys(s::MinMaxInterpStack) = keys(specs(s)) 

# Base.getindex(s::MinMaxInterpStack, key::Symbol, i1::StandardIndices, I::StandardIndices...) =
#     if key in keys(s.minmaxinterpmodes)
#         arrays = 
#         arrays = map(a -> readwindowed(a, i1, I...), arrays)
#         InterpArray(arrays, int[key], s.frac)
#     else
#         interparray(stack, key, i1, I...)
#     end

function minmaxseries(series::GeoSeries, dates, interpolators, specs::NamedTuple)
    # Convert series to mutable ReadOnceStack that will become GeoStacks 
    # the first time they are accessed. Interpolated sliced will share 
    # CachedStack so they are only loaded once.
    origdates = index(series, Ti)
    cseries = CachedStack.(series)
    interpseries = Vector(undef, length(dates))
    for (i, t) in enumerate(dates)
        # Find the t in the serie Ti index 
        x = searchsortedlast(index(cseries, Ti), t)
        # We need 3 stacks, as we need the max temp
        # from the previous day for times before tmin.
        stacks = if x <= firstindex(cseries)
            error("Include dates at least 1 period before the first in the required dates")
        elseif x == lastindex(cseries)
            error("Include dates at least 1 period beyond the last in the required dates")
        else
            [cseries[x-1], cseries[x], cseries[x+1]]
        end
        stacks = OffsetArray(stacks, 0:2)
        frac = calcfrac(origdates[x], origdates[x] + step(origdates), t)
        minmaxfracs = map(spec -> MinMaxFracs(spec, step(dates), frac), specs)
        interpseries[i] = MinMaxInterpStack(stacks, interpolators, frac, minmaxfracs)
    end
    GeoSeries(interpseries, Ti(dates); childtype=InterpStack)
end

function meandayminmaxseries(series::GeoSeries, dates, interpolators, specs::NamedTuple)
    # Convert series to mutable CachedStack that will become GeoStacks 
    # the first time they are accessed. Interpolated sliced will share 
    # CachedStack so they are only loaded once.
    origdates = index(series, Ti)
    cseries = CachedStack.(series)
    interpseries = Vector(undef, length(dates))
    for (i, t) in enumerate(dates)
        # Find the t in the series Ti index 
        x = searchsortedlast(index(cseries, Ti), t)
        # We need 3 stacks, as we need the max temp
        # from the previous day for times before tmin.
        stacks = if x < firstindex(cseries) || x > lastindex(cseries)
            error("Date $x is outside the series")
        else
            # We just use the same average day three times
            [cseries[x], cseries[x], cseries[x]]
        end
        stacks = OffsetArray(stacks, 0:2)
        frac = calcfrac(origdates[x], origdates[x] + step(origdates), t)
        minmaxfracs = map(spec -> MinMaxFracs(spec, step(dates), frac), specs)
        interpseries[i] = MinMaxInterpStack(stacks, interpolators, frac, minmaxfracs)
    end
    dim = if dates isa Dimension
        dates
    else
        Ti(dates)
    end
    GeoSeries(interpseries, dim; childtype=InterpStack)
end

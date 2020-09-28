struct MinMaxSpec{T,IM}
    times::T
    interpmode::IM
end

interpmode(mm::MinMaxSpec) = mm.interpmode

function windowfromtimes(spec::MinMaxSpec, step, tfrac) 
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
    inds, frac = windowfromtimes(spec, step, trfrac)
    imode = spec.interpmode
    MinMaxFracs{typeof(inds),typeof(frac),typeof(imode)}(inds, frac, imode)
end

interpmode(mm::MinMaxFracs) = mm.interpmode


"""
    MimMaxInterpStack

Provides interpolation for a value (usually temperature) that usually
has minimum and maximum values provided in separate layers of the
stack. Other layers in the stack are interpolated, or not,
as for `InterpStack`.
"""
struct MinMaxInterpStack{S,I,MMS<:NamedTuple{<:Any,Tuple{Vararg{<:MinMaxFracs}}}} <: AbstractGeoStack{S}
    stacks::S
    interpmodes::I
    frac::Float64
    specs::MMS
end
MinMaxInterpStack(stacks, interpmodes, tfrac, specs::NamedTuple{<:Any,<:Tuple{Vararg{<:MinMaxSpec}}}, timestep) = begin
    minmaxfracs = map(spec -> MinMaxFracs(spec, timestep, tfrac), specs)
    interpseries[i] = MinMaxInterpStack(stacks, interpmodes, frac, minmaxfracs)
end

interpmodes(s::MinMaxInterpStack) = s.interpmodes
frac(s::MinMaxInterpStack) = s.frac
specs(s::MinMaxInterpStack) = s.specs

Base.getindex(s::MinMaxInterpStack, key::Symbol) =
    if key in keys(s.specs)
        spec = s.specs[key]
        intmode = interpmode(spec)
        arrays = map((k, i) -> s.stacks[i][k], keys(spec.indices), spec.indices)
        imode = interpmode(spec)
        valweights = Interpolations.value_weights(i.degree, frac)
        weights = Interpolations.WeightedIndex(1:2, valweights)
        InterpArray(arrays, imode, frac, interpolate(SVector(arrays), imode), weights)
    else
        interparray(stack, key)
    end

# Base.getindex(s::MinMaxInterpStack, key::Symbol, i1::StandardIndices, I::StandardIndices...) =
#     if key in keys(s.minmaxinterpmodes)
#         arrays = 
#         arrays = map(a -> readwindowed(a, i1, I...), arrays)
#         InterpArray(arrays, int[key], s.frac)
#     else
#         interparray(stack, key, i1, I...)
#     end

function minmaxseries(series::GeoSeries, newtimeseries, interpolators, spec::NamedTuple)
    # Convert series to mutable ReadOnceStack that will become GeoStacks 
    # the first time they are accessed. Interpolated sliced will share 
    # ReadOnceStack so they are only loaded once.
    origtimeseries = index(series, Ti)
    ro_series = ReadOnceStack.(series)
    T = Union{eltype(ro_series),Missing}
    interpseries = Vector(undef, length(newtimeseries))
    for (i, t) in enumerate(newtimeseries)
        # Find the t in the serie Ti index 
        x = searchsortedlast(index(ro_series, Ti), t)
        # We need 3 stack, as we need the max temp
        # from the previous day for times before tmin.
        stacks = if x < firstindex(ro_series)
            T[missing, missing, ro_series[x]]
        elseif x == firstindex(ro_series)
            T[missing, ro_series[x], ro_series[x+1]]
        elseif x == lastindex(ro_series)
            T[ro_series[end-1], ro_series[end], missing]
        else
            T[ro_series[x-1], ro_series[x], ro_series[x+1]]
        end
        stacks = OffsetArray(stacks, 0:2)
        frac = calcfrac(origtimeseries[x], origtimeseries[x], t)
        minmaxfracs = map(MinMaxFracs(spec, step(newtimeseries), frac), specs)
        interpseries[i] = MinMaxInterpStack(stacks, interpolators, frac, minmaxfracs)
    end
    GeoSeries(interpseries, Ti(newtimeseries); childtype=InterpStack)
end



# Interpolation methods for cyclic min/max values like temperature

"""
    Cosine()

Similar to a `Linear` interpolator, but using `sin` to calculate
the weights, as if it is cyclic between two values.
"""
struct Cosine <: Interpolations.Degree{1} end

Interpolations.positions(::Cosine, ax::AbstractUnitRange{<:Integer}, x) =
    Interpolations.positions(Linear(), ax, x)

Interpolations.value_weights(::Cosine, δx) = begin
    # calculate the fraction from a sin
    cosfrac = -cos(π * δx)
    # normalise between 0 and 1
    normed = (cosfrac + 1) / 2
    (1-normed, normed)
end
Interpolations.gradient_weights(::Cosine, δx) = (-oneunit(δx), oneunit(δx))
Interpolations.hessian_weights(::Cosine, δx) = (zero(δx), zero(δx))
Interpolations.padded_axis(ax::AbstractUnitRange, ::BSpline{Cosine}) = ax

"""
    HypTan{S}()

Similar to a `Linear` interpolator, but using hyperbolic tangent to 
calculate the weights. It has a multiplier `X` that controls how long the
curve is near each index. This is automatically scaled so the max and min 
values are always 0 and 1 for all `X`. At the limits, HypTan is essentially 
`Linear` and `NoInterp`. 

`HypTan{0}` will error, instead uses `Linear` - which is identical.
"""
struct HypTan{X,S} <: Interpolations.Degree{1} end
HypTan() = HypTan{2}()
HypTan{X}() where X = HypTan{X,1 / tanh(π * X / 2)}()
HypTan{0}() = Linear()

Interpolations.positions(::HypTan, ax::AbstractUnitRange{<:Integer}, x) =
    Interpolations.positions(Linear(), ax, x)

Interpolations.value_weights(ht::HypTan{X,S}, δx) where {X,S} = begin
    # calculate the fraction from a sin
    hypfrac = tanh(X * π * (δx - 0.5)) * S
    # normalise between 0 and 1
    normed = (hypfrac + 1) / 2
    (1-normed, normed)
end
Interpolations.gradient_weights(::HypTan, δx) = (-oneunit(δx), oneunit(δx))
Interpolations.hessian_weights(::HypTan, δx) = (zero(δx), zero(δx))
Interpolations.padded_axis(ax::AbstractUnitRange, ::BSpline{HypTan}) = ax

"""
    MinMaxIterpolator(times, interptype)

Specifies the time that correspond to the daily minimum and maximum.

# Arguments

- `times`: a length 2 `NamedTuple` where the keys correspond to the layer names for the
    minumum and maximum values, and the values specify the offset in hours at which they occur.
- `interptype`: An Interpolations.jl `InterpType`. This package provides [`Cosine`](@ref)
    and [`Hyptan`](@ref).
"""
struct MinMaxInterpolator{T,IM}
    times::T
    interptype::IM
end

interptype(mm::MinMaxInterpolator) = mm.interptype
Base.keys(mm::MinMaxInterpolator) = Base.keys(mm.times)

# Used internally, replacing MinMaxInterpolator with
# interpolators for specific interpolation fractions.
struct MinMaxFrac{I,F,IM}
    indices::I
    frac::F
    interptype::IM
end
function MinMaxFrac(interpolator::MinMaxInterpolator, step, tfrac::AbstractFloat)
    key_indices, frac = chooseinterparrays(interpolator, step, tfrac)
    itype = interpolator.interptype
    return MinMaxFrac{typeof(key_indices),typeof(frac),typeof(itype)}(key_indices, frac, itype)
end

function chooseinterparrays(interpolator, step, tfrac)
    # Get layers, indices and interpolation fraction for this specific
    # fraction of the timestep
    minfrac, maxfrac = map(t -> Second(t) / Second(step), interpolator.times)
    # The last section, larger than the max time
    if tfrac >= maxfrac
        stackkeys = reverse(keys(interpolator.times))
        window = maxfrac, minfrac + 1
        inds = 1, 2
    # The first section
    elseif tfrac < minfrac
        stackkeys = reverse(keys(interpolator.times))
        window = maxfrac - 1, minfrac
        inds = 0, 1
    # The middle section
    else
        stackkeys = keys(interpolator.times)
        window = minfrac, maxfrac
        inds = 1, 1
    end
    frac = -(window[1] - tfrac) / (window[2] - window[1])
    key_indices = NamedTuple{stackkeys}(inds)
    return key_indices, frac
end

frac(mm::MinMaxFrac) = mm.frac
indices(mm::MinMaxFrac) = mm.indices
interptype(mm::MinMaxFrac) = mm.interptype

Base.keys(mm::MinMaxInterpolator) = Base.keys(mm.times)

function minmaxarray(stacks::OffsetArray, mm_fracs::NamedTuple, key::Symbol)
    if key in keys(mm_fracs)
        # Min/max interpolation key
        mm_frac = mm_fracs[key]
        frac = mm_frac.frac
        itype = interptype(mm_frac)
        layers = map((k, i) -> stacks[i][k], keys(mm_frac.indices), mm_frac.indices)
        if frac < 1
            frac += 1
        end
        A1 = first(layers)
        data = map(parent, layers)
        iA = InterpArray(data, itype, frac)
        return GeoArray(iA, dims(A1); name=key, missingval=missingval(A1))
    else
        # Other key
        return interparray(stacks, interptypes, frac, key)
    end
end

"""
    minmaxstack(stacks, interptypes, tfrac, mm_fracs)
    minmaxstack(stacks, interptypes, tfrac, mm_interpolators, timestep)

An `AbstractGeoStack` that Provides interpolation for a value (often temperature)
that usually has minimum and maximum values provided in separate layers.
Other layers in the stack are interpolated, or not, as for [`InterpStack`](@ref).

# Arguments

- `stack`: a `Tuple` of `AbstractGeoStack`, usually of length 3.
- `interptypes`: a `NamedTuple` connecting stack layer name with an Interpolations.jl `InterpType`.
    These are used for other layers used with the min/max layers.
- `tfrac`: the fraction of the outer timestep this slice is taken at
- `mm_fracs`: a `NamedTuple` connecting new layer names (e.g. `:temp`) to a `MinMaxFrac` using other stack layers.
- `mm_interpolators`: a `NamedTuple` connecting new layer names (e.g. `:temp`) to a `MinMaxFrac` using other stack layers.
- `timestep`: the size of the inner timestep, e.g `Hour(1)`
"""
function minmaxstack(
    stacks, interptypes, frac::Float64, mm_interpolators::NamedTuple, step,
)
    # Combine of keys defined in mm_interpolators and the other
    # layers in the stack that are not use in min/max interpolation.
    minmaxkeys = keys(mm_interpolators)
    minmax_sourcekeys = Tuple(Iterators.flatten(map(keys, mm_interpolators)))
    stackkeys = keys(first(stacks))
    combinedkeys = foldl(stackkeys; init=minmaxkeys) do acc, key
        key in minmax_sourcekeys ? acc : (acc..., key)
    end
    # Get specific fractional interpolation objects for this tfrac/step combination
    # Usually there will only be one interpolator - but it needs a key anyway so we
    # use a NamedTuple and `map`
    mm_fracs = map(mm_interpolators) do interp
        MinMaxFrac(interp, step, frac)
    end
    # Generate InterpArray layers for each key in the stack
    layers = map(combinedkeys) do key
        if key in keys(mm_fracs)
            minmaxarray(stacks, mm_fracs, key)
        else
            interparray(stacks, interptypes, frac + 1, key)
        end
    end
    return GeoStack(layers)
end

"""
    minmaxseries(series::GeoSeries, dates; interptypes, mm_interpolators::NamedTuple)

Generate a `GeoSeries` of `MinMaxStack` from `series` for each date/time in `dates`.
Min/max interpolators are defined in the `mm_interpolators`, a `NamedTuple` of

`dates` must define `step`. An `AbstractRange` or a DimensionalData.lj `Dimension` can also be used.

# Arguments

- `series`: a `GeoSeries` of `GeoStack`, with a `Regular` `Sampled` `Ti` dimension.
- `dates`: `AbstractVector` of dates
- `interptypes`: a `NamedTuple` connecting stack layer name with an Interpolations.jl `InterpType`.
    These are used for other layers used with the min/max layers.
- `mm_interpolators`: a `NamedTuple` connecting new layer names (e.g. `:temp`) to a `MinMaxFrac` using other stack layers.
"""
function minmaxseries(series::GeoSeries, dates;
    step=step(dates), interptypes=NamedTuple(), mm_interpolators::NamedTuple
)
    hasdim(series, Ti) || throw(ArgumentError("series must have a Ti dimension"))
    origdates = dims(series, Ti)
    interpseries = map(dates) do t
        # Find the t in the serie Ti index
        x = searchsortedlast(index(series, Ti), t)
        # We need 3 stacks, as we need the max temp
        # from the previous day for times before tmin.
        stacks = if x <= firstindex(series)
            error("Include dates at least 1 period before the first in the required dates")
        elseif x == lastindex(series)
            error("Include dates at least 1 period beyond the last in the required dates")
        else
            [series[x-1], series[x], series[x+1]]
        end
        stacks = OffsetArray(stacks, 0:2)
        frac = calcfrac(origdates[x], origdates[x] + Base.step(origdates), t)
        minmaxstack(stacks, interptypes, frac, mm_interpolators, Base.step(origdates))
    end
    dim = timedim(dates, step)
    return GeoSeries(interpseries, dim)
end

"""
    meanday_minmaxseries(series::GeoSeries, dates, interptypes, mm_interpolators::NamedTuple) => GeoSeries

Generate a `GeoSeries` of `MinMixStack` for one day in each period in `series`, with slices
for every `step` in the day.

# Arguments

- `series`: a `GeoSeries` of `GeoStack`
- `dates`: `AbstractVector` of dates
- `step`: the step size of the intervals in the resulting `GeoSeries`.
- `interptypes`: a `NamedTuple` connecting stack layer name with an Interpolations.jl `InterpType`.
    These are used for other layers used with the min/max layers.
- `mm_interpolators`: a `NamedTuple` connecting new layer names (e.g. `:temp`) to a `MinMaxFrac` using other stack layers.
"""
function meanday_minmaxseries(series::GeoSeries, dates;
    step=step(dates), interptypes=NamedTuple(), mm_interpolators::NamedTuple
)
    hasdim(series, Ti) || throw(ArgumentError("series must have a Ti dimension"))
    origdates = dims(series, Ti)
    interpseries = map(dates) do t
        # Find the t in the series Ti index
        x = searchsortedlast(index(series, Ti), t)
        # We need 3 stacks, as we need the max temp
        # from the previous day for times before tmin.
        if t < DD.bounds(series, Ti)[1] || t > DD.bounds(series, Ti)[2]
            error("Date $x is outside the series")
        end
        # We just use the same average day three times
        stacks = OffsetArray([series[x], series[x], series[x]], 0:2)
        daydate = DateTime(yearmonthday(t)...)
        dayfrac = calcfrac(daydate, daydate + Day(1), t)
        minmaxstack(stacks, interptypes, dayfrac, mm_interpolators, Day(1))
    end
    dim = timedim(dates, step)
    return GeoSeries(interpseries, dim)
end

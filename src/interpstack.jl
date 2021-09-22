"""
    InterpArray(arrays::Tuple, fraction::Tuple)

InterplatedArray holds multiple arrays the will be interpolated when indexed.
The fractions to include from each arrays are captured in the `fractions` argument.
You may use as many arrays as required.
"""
struct InterpArray{T,N,A<:AbstractVector{<:AbstractArray{<:Any,N}},B,F,W} <: DA.AbstractDiskArray{T,N}
    arrays::A
    interptype::B
    frac::F
    weights::W
end
InterpArray(arrays::Tuple, interptype, frac) = InterpArray([arrays...], interptype, frac)
function InterpArray(arrays::AbstractVector, interptype, frac)
    A1 = first(arrays)
    valweights = Interpolations.value_weights(interptype.degree, frac)
    weights = Interpolations.WeightedIndex(axes(arrays, 1), valweights)
    arrays = [arrays...]
    A1 = first(arrays)
    # Promote to Float64 as we will interpolate the value
    T = eltype(A1) isa AbstractFloat ? eltype(A1) : Float64
    N = ndims(A1)
    InterpArray{T,N,typeof.((arrays, interptype, frac, weights))...}(
        arrays, interptype, frac, weights
    )
end

arrays(A::InterpArray) = A.arrays
interptype(A::InterpArray) = A.interptype
frac(A::InterpArray) = A.frac
weights(A::InterpArray) = A.weights

DA.haschunks(A::InterpArray) = DA.haschunks(first(arrays(A)))
DA.eachchunk(A::InterpArray) = DA.eachchunk(first(arrays(A)))
DA.readblock!(A::InterpArray, aout, r::AbstractUnitRange...) = aut .= A[r...]

# Base methods

# Copy with a broadcast so interpolation will work on a GPU?
Base.copy!(dst::AbstractArray, src::InterpArray) = dst .= src

Base.size(A::InterpArray) = Base.size(first(A.arrays))
Base.size(A::InterpArray, dims::Int) = Base.size(first(A.arrays), dims)
function Base.getindex(A::InterpArray, i1::StandardIndices, I::StandardIndices...)
    data = [map(a -> a[i1, I...], arrays(A))...]
    int = interpolate(data, interptype(A))
    return collect(int(frac(A)))
end
function Base.getindex(A::InterpArray, i1::Int, I::Int...)
    Vector(map(a -> a[i1, I...], arrays(A)))[A.weights]
end

function interparray(stacks, interptypes, frac, key::Symbol)
    layers = map(p -> p[key], stacks)
    A1 = first(layers)
    if haskey(interptypes, key)
        A = InterpArray(map(parent, layers), interptypes[key], frac)
        return GeoArray(A, dims(A1); name=key, missingval=missingval(A1))
    else
        return centerval(layers)
    end
end

function interpstack(stacks, interptypes, frac)
    stackkeys = keys(first(stacks))
    data = map(stackkeys) do key
        parent(interparray(stacks, interptypes, frac + 1, key))
    end |> NamedTuple{stackkeys}
    return GeoStack(first(stacks); data=data)
end

centerval(vals) = vals[length(vals) ÷ 2 + 1]

# Series loading code

nlayers(bspline::BSpline) = nlayers(bspline.degree)
nlayers(::NoInterp) = 1
nlayers(::Interpolations.Degree{X}) where X = X

function interpseries(series::GeoSeries, dates; interptypes, step=step(dates))
    # Convert series to mutable CachedStack that will become GeoStacks 
    # the first time they are accessed. Interpolated sliced will share 
    # CachedStack so they are only loaded once.
    origtimeseries = index(series, Ti)
    T = Union{eltype(series),Missing}
    interpseries = map(dates) do t
        # Find the t in the serie Ti index 
        x = searchsortedlast(index(series, Ti), t)
        # For now we only allow length degree 2 interptypes,
        # using 2 stacks.
        stacks = if x <= firstindex(series)
            error("Include dates at least 1 period before the first in the required dates")
        elseif x >= lastindex(series)
            error("Inlcude dates at least 1 period beyond the last in the required dates")
        else
            T[series[x], series[x+1]]
        end
        frac = calcfrac(origtimeseries[x], origtimeseries[x] + Base.step(origtimeseries), t)
        interpstack(stacks, interptypes, frac)
    end
    return GeoSeries(interpseries, timedim(dates, step))
end

timedim(dates::Dimension, step) = dates
function timedim(dates::AbstractRange, step)
    Ti(dates; mode=Sampled(order=Ordered(), span=Regular(step), sampling=Intervals(Start())))
end
function timedim(dates::AbstractArray, step)
    span = Explicit(permutedims(hcat(dates, dates .+ step)))
    return Ti(dates; mode=Sampled(order=Ordered(), span=span, sampling=Intervals(Start())))
end

calcfrac(a, b, x) = (x - a) / (b - a)
calcfrac(a::DateTime, b::DateTime, x::DateTime) = (x - a).value / (b - a).value

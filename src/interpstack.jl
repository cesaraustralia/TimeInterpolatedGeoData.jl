"""
    InterpArray(arrays::Tuple, fraction::Tuple)

InterplatedArray holds multiple arrays the will be interpolated
when indexed. The fractions to include from each arrays are
captured in the `fractions` argument. You may use as many arrays
as required.
"""
struct InterpArray{T,N,D,A<:AbstractArray{<:Union{Missing,AbstractGeoArray{T,N,D}}},B,F,I,W} <: AbstractGeoArray{T,N,D,A}
    arrays::A
    interpmode::B
    frac::F
    interpolator::I
    weights::W
end
InterpArray(arrays, interpmode, frac) = begin
    valweights = Interpolations.value_weights(interpmode.degree, frac)
    weights = Interpolations.WeightedIndex(1:2, valweights)
    arrays = [arrays...]
    InterpArray(arrays, interpmode, frac, interpolate(arrays, interpmode), weights)
end

arrays(A::InterpArray) = A.arrays
interpmode(A::InterpArray) = A.interpmode
frac(A::InterpArray) = A.frac
interpolator(A::InterpArray) = A.interpolator
weights(A::InterpArray) = A.weights

GeoData.data(A::InterpArray) = interpolator(A)[frac(A)]
GeoData.dims(A::InterpArray) = dims(first(arrays(A)))
GeoData.refdims(A::InterpArray) = refdims(first(arrays(A))) 
GeoData.name(A::InterpArray) = GeoData.name(first(arrays(A)))
GeoData.missingval(A::InterpArray) = missingval(first(arrays(A)))
GeoData.metadata(A::InterpArray) = metadata(first(arrays(A)))

# Base methods

# Copy with a broadcast so interpolation will work on a GPU?
Base.copy!(dst::AbstractArray, src::InterpArray) = dst .= src

Base.size(s::InterpArray, args...) = size(first(arrays(s)), args...)
Base.getindex(A::InterpArray, I::StandardIndices...) = begin
    As = [map(a -> view(a, I...), arrays(A))...]
    int = interpolate(As, interpmode(A))
    int(frac(A))[I...]
end
Base.getindex(A::InterpArray, i1::Int, I::Int...) = 
    Vector(map(a -> a[i1, I...], arrays(A)))[A.weights]

interparray(s, key::Symbol) = begin
    arrays = map(p -> p[key], stacks(s))
    intmodes = interpmodes(s)
    if haskey(intmodes, key)
        InterpArray(arrays, intmodes[key], s.frac)
    else
        centerval(arrays)
    end
end
interparray(s, key::Symbol, i1::StandardIndices, I::StandardIndices...) = begin
    arrays = map(p -> p[key], stacks(s)) 
    intmodes = interpmodes(s) 
    if haskey(intmodes, key)
        arrays = map(a -> readwindowed(a, i1, I...), arrays)
        InterpArray(arrays, intmodes[key], s.frac)
    else
        readwindowed(centerval(arrays), i1, I...)
    end
end

struct InterpStack{P,I<:NamedTuple} <: AbstractGeoStack{P}
    stacks::P
    interpmodes::I
    frac::Float64
end
stacks(s::InterpStack) = s.stacks
interpmodes(s::InterpStack) = s.interpmodes
frac(s::InterpStack) = s.frac

centerval(vals) = vals[(length(vals) + 1) ÷ 2 + 1]
centerval(s::InterpStack) = centerval(stacks(s))

# GeoData methods

GeoData.refdims(s::InterpStack) = refims(centerval(s))
GeoData.metadata(s::InterpStack) = metadata(centerval(s))
GeoData.window(s::InterpStack) = window(centerval(s))
GeoData.childtype(s::InterpStack) = childtype(centerval(s))

# Base methods

Base.keys(s::InterpStack) = keys(centerval(s))
Base.keys(s::CachedStack) = keys(centerval(s))
Base.values(s::CachedStack) = values(centerval(s))
Base.names(s::CachedStack) = names(centerval(s))
Base.length(s::AbstractGeoStack) = length(centerval(s))

Base.getindex(s::InterpStack, key::Symbol) = interparray(s, key)
Base.getindex(s::InterpStack, key::Symbol, i1::StandardIndices, I::StandardIndices...) =
    interparray(s, key, i1, I...)

clear!(s::InterpStack) = (map(clear!, s.stacks); nothing)



# Series loading code

nlayers(bspline::BSpline) = nlayers(bspline.degree)
nlayers(::NoInterp) = 1
nlayers(::Interpolations.Degree{X}) where X = X

function interpseries(series::GeoSeries, newtimeseries, interpolators)
    # Convert series to mutable CachedStack that will become GeoStacks 
    # the first time they are accessed. Interpolated sliced will share 
    # CachedStack so they are only loaded once.
    origtimeseries = index(series, Ti)
    ro_series = CachedStack.(series)
    T = Union{eltype(ro_series),Missing}
    interpseries = Vector(undef, length(newtimeseries))
    for (i, t) in enumerate(newtimeseries)
        # Find the t in the serie Ti index 
        x = searchsortedlast(index(ro_series, Ti), t)
        # For now we only allow length degree 2 interpolators,
        # using 2 stacks.
        stacks = if x <= firstindex(ro_series)
            error("Include dates at least 1 period before the first in the required dates")
        elseif x >= lastindex(ro_series)
            error("Inlcude dates at least 1 period beyond the last in the required dates")
        else
            T[ro_series[x], ro_series[x+1]]
        end
        frac = calcfrac(origtimeseries[x], origtimeseries[x], t)
        interpseries[i] = InterpStack(stacks, interpolators, frac)
    end
    GeoSeries(interpseries, Ti(newtimeseries); childtype=InterpStack)
end

calcfrac(a, b, x) = (x - a) / (b - a)
calcfrac(a::DateTime, b::DateTime, x::DateTime) = (x - a).value / (b - a).value

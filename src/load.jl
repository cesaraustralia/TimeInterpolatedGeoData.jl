
# Define single array and time-series loaders for WorldClim and AWAP
# These will move to another package at some stage, it's just not clear 
# where yet.

# WorldClim Weather

geoarray(T::Type{<:WorldClim}, args...; kwargs...) = 
    GDALarray(rasterpath(T, args...); kwargs...)

stack(T::Type{WorldClim{Weather}}; date, layers=rasterlayers(T), kwargs...) = 
    GDALstack(map(l -> rasterpath(T, l, date), layers); keys=layers, kwargs...)

series(T::Type{WorldClim{Weather}}; dates, layers=rasterlayers(T), window=(), kwargs...) = begin
    dates = first(dates):Month(1):last(dates)
    timedim = Ti(dates)
    stacks = [stack(T; date=d, layers=layers, window=window) for d in dates]
    GeoData.GeoSeries(stacks, timedim; kwargs...)
end


# ALWB

function geoarray(T::Type{<:ALWB}, layer::Type, date::TimeType; kwargs...)
    key = Symbol(RasterDataSources._pathsegment(layer))
    NCDarray(rasterpath(T, layer, date), key; kwargs...)
end

stack(T::Type{<:ALWB}; date, layers=rasterlayers(T), kwargs...) = 
    GeoStack(map(l -> geoarray(T, l, date), layers)...; kwargs...)

series(T::Type{<:ALWB}; dates, layers=rasterlayers(T), window=(), kwargs...) = begin
    dates = first(dates):Month(1):last(dates)
    timedim = Ti(dates)
    stacks = [stack(T; date=d, layers=layers, window=window) for d in dates]
    GeoData.GeoSeries(stacks, timedim; kwargs...)
end


# AWAP

# AWAP asciii files don't have crs and BOM don't specify what it is besides
# being evenly spaced lat/lon. So we assume it's EPSG(4326) and set it manually

geoarray(T::Type{AWAP}, args...; kwargs...) = 
    GDALarray(rasterpath(T, args...); crs=EPSG(4326), kwargs...)

stack(T::Type{AWAP}; date, layers=rasterlayers(T), kwargs...) = 
    GDALstack(map(l -> rasterpath(T, l, date), layers); childkwargs=(crs=EPSG(4326),), 
              keys=map(layerkey, layers), kwargs...)

series(T::Type{AWAP}; dates, layers=rasterlayers(T), window=(), kwargs...) = begin
    dates = first(dates):Day(1):last(dates)
    timedim = Ti(dates)
    stacks = [stack(T; date=d, layers=layers, window=window) for d in dates]
    GeoData.GeoSeries(stacks, timedim; kwargs...)
end

# TODO these keys are pretty bad
layerkey(t) = name(typeof(t))
layerkey(::Type{<:Solar}) = :solar
layerkey(::Type{<:Rainfall}) = :rainfall
layerkey(::Type{VapourPressure{H09}}) = :vpr09
layerkey(::Type{VapourPressure{H15}}) = :vpr15
layerkey(::Type{Temperature{MinAve}}) = :tmin
layerkey(::Type{Temperature{MaxAve}}) = :tmax


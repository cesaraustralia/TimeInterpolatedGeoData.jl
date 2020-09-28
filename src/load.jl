
# Define single array and time-series loaders for WorldClim and AWAP
# These will move to another package at some stage, it's just not clear 
# where yet.

# Worldclim uses symbol keys for the layer
geoarray(T::Type{<:WorldClim}, args...; kwargs...) = 
    GDALarray(rasterpath(T, args...); kwargs...)

series(T::Type{WorldClim{Weather}}, layers=WEATHER_LAYERS; window=(), kwargs...) = begin
    namedlayers = NamedTuple{layers}(layers)
    layerfiles = map(namedlayers) do layer
        readdir(rasterpath(T, layer); join=true)
    end
    all(map(l -> length(l) == length(first(layerfiles)), layerfiles)) || error("Wrong number of files")
    firstdate = SSL._filename2date(T, first(layerfiles[1]))
    lastdate = SSL._filename2date(T, last(layerfiles[1]))
    timedim = Ti(firstdate:Month(1):lastdate)
    stacks = [DiskStack(map(fns -> fns[i], layerfiles); window=window, childtype=GDALarray) for i in eachindex(first(layerfiles))]
    GeoData.GeoSeries(stacks, timedim; kwargs...)
end

# AWAAP uses types to specify the layer. 
# I'm not sure which makes the most sense
geoarray(T::Type{AWAP}, t, date; kwargs...) = 
    GDALarray(rasterpath(T, t, date); kwargs...)

series(T::Type{AWAP}, layers=SSL.AWAP_LAYERS; kwargs...) = begin
    T = AWAP
    keys = map(layerkey, layers)
    namedlayers = NamedTuple{keys}(layers)
    layerfiles = map(namedlayers) do layer
        readdir(rasterpath(AWAP, layer); join=true)
    end
    all(map(l -> length(l) == length(first(layerfiles)), layerfiles)) || error("Wrong number of files")
    firstdate = splitext(basename(first(layerfiles[1])))[1] 
    lastdate = splitext(basename(last(layerfiles[1])))[1]
    timedim = Ti(SSL._string2date(T, firstdate):Day(1):SSL._string2date(T, lastdate))
    stacks = [DiskStack(map(fns -> fns[i], layerfiles); childtype=GDALarray) for i in eachindex(first(layerfiles))]
    GeoData.GeoSeries(stacks, timedim; kwargs...)
end

layerkey(t) = name(typeof(t))
layerkey(::Type{<:Solar}) = :solar
layerkey(::Type{<:Rainfall}) = :rainfall
layerkey(::Type{VapourPressure{H09}}) = :vprph09
layerkey(::Type{VapourPressure{H15}}) = :vprph15
layerkey(::Type{Temperature{MinAve}}) = :temp_minave
layerkey(::Type{Temperature{MaxAve}}) = :temp_maxave


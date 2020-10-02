using TimeInterpolatedGeoData, SimpleSDMDataSources, Test, Dates, GeoData

using SimpleSDMDataSources: Values, SoilMoisture, Upper, Lower

@testset "load WorldClim Weather" begin
    # Weather
    rasterpath(WorldClim{Weather}, :prec, Date(2001, 1))
    download_raster(WorldClim{Weather}, :prec,; dates=Date(2001):Month(1):Date(2001, 12))

    # Weather time-series
    dates = (Date(2001), Date(2002))
    ser = series(WorldClim{Weather}; dates=dates, window=(Lat(600:1900), Lon(600:1900)))
    ser[Date(2001, 1)][:tmax]
    # Select Australia, using regular lat/lon selectors
    A = geoarray(WorldClim{Weather}, :prec, Date(2001, 05); usercrs=EPSG(4326))
    A[Lat(Between(-10, -45)), Lon(Between(110, 160))]
end

@testset "load ALWB" begin
    dates = DateTime(2019, 09, 19), DateTime(2019, 11, 19)
    layers = (SoilMoisture{Upper}, SoilMoisture{Lower},) 
    download_raster(ALWB{Values,Day}, layers; dates=dates)
    st = stack(ALWB{Values,Day}; layers=layers, date=DateTime(2019, 10, 19))
    A = geoarray(ALWB{Values,Day}, SoilMoisture{Lower}, DateTime(2019, 10, 19))
end

@testset "load AWAP" begin
    dates = DateTime(2019, 09, 19), DateTime(2019, 11, 19)
    download_raster(AWAP; dates=dates)
    A = geoarray(AWAP, Rainfall, DateTime(2019, 10, 19))
    A = geoarray(AWAP, Temperature{MaxAve}, DateTime(2019, 10, 19))
    series(AWAP; dates=dates)[DateTime(2019, 10, 1)][:solar]
end

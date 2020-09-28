using TimeInterpolatedGeoData

# Load and plot

# Weather
rasterpath(WorldClim{Weather}, :prec, Date(2001, 1))
SimpleSDMDataSources._filename2date(WorldClim, "wc2.1_2.5m_tmax_2001-02.tif")
download_raster(WorldClim{Weather}; layers=(:prec,), dates=Date(2001):Month(1):Date(2001, 12))

geoarray(WorldClim{Weather}, :prec, Date(2001)) |> plot


# BioClim
rasterpath(WorldClim{BioClim}, 1, 10.0)
download_raster(WorldClim{BioClim}; layer=1, resolution=10.0)

geoarray(WorldClim{BioClim}, 1, 10.0) |> plot


# Weather time-series
ser = series(WorldClim{Weather}; window=(Lat(600:1900), Lon(600:1900)))
ser[Date(2001, 1)][:tmax] |> plot
# Plot Australia, using regular lat/lon selectors
A = geoarray(WorldClim{Weather}, :prec, Date(2001, 05); usercrs=EPSG(4326))
A[Lat(Between(-10, -45)), Lon(Between(110, 160))] |> plot


# AWAP
dates = DateTime(2019, 09, 19), DateTime(2019, 11, 19)
download_raster(AWAP; dates=dates)
A = geoarray(AWAP, Rainfall, DateTime(2019, 10, 19))
plot(A)
A = geoarray(AWAP, Temperature{MaxAve}, DateTime(2019, 10, 19))
plot(A)
series(AWAP)
[DateTime(2019, 10, 1)][:solar]
|> plot

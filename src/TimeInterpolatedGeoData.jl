module TimeInterpolatedGeoData

using GeoData, 
      Interpolations, 
      Dates, 
      StaticArrays, 
      OffsetArrays, 
      RasterDataSources,
      ArchGDAL, 
      NCDatasets

using RasterDataSources: Solar, Rainfall, VapourPressure, Temperature, 
      H09, H15, MinAve, MaxAve, layers

using GeoData: StandardIndices, readwindowed

export CachedStack 

export InterpStack, InterpArray

export MinMaxInterpStack, MinMaxSpec

export Cosine, HypTan

export geoarray, stack, series, interpseries, minmaxseries, meandayminmaxseries 

include("cachedstack.jl")
include("interpstack.jl")
include("minmaxstack.jl")
include("interpolations.jl")
include("load.jl")

end # module

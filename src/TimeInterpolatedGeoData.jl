module TimeInterpolatedGeoData

using GeoData, 
      Interpolations, 
      Dates, 
      StaticArrays, 
      OffsetArrays, 
      RasterDataSources,
      ArchGDAL, 
      NCDatasets

using GeoData: StandardIndices, readwindowed

export CachedStack 

export InterpStack, InterpArray

export MinMaxInterpStack, MinMaxSpec

export Cosine, HypTan

export interpseries, minmaxseries, meandayminmaxseries

include("cachedstack.jl")
include("interpstack.jl")
include("minmaxstack.jl")
include("interpolations.jl")

end # module

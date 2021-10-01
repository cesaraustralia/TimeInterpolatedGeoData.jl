module TimeInterpolatedGeoData

using GeoData, 
      Interpolations, 
      Dates, 
      StaticArrays, 
      OffsetArrays, 
      RasterDataSources

using GeoData: StandardIndices

export CachedStack, InterpArray, MinMaxInterpolator

export Cosine, HypTan

export interpstack, minmaxinterpstack 

export interpseries, minmaxseries, meanday_minmaxseries

const GD = GeoData
const DD = GeoData.DimensionalData
const DA = GeoData.DiskArrays

include("cachedstack.jl")
include("interpstack.jl")
include("minmaxstack.jl")
include("interpolations.jl")

end # module

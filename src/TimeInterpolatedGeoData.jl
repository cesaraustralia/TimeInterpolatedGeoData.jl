module TimeInterpolatedGeoData

using GeoData, 
      Interpolations, 
      Dates, 
      StaticArrays, 
      OffsetArrays, 
      SimpleSDMDataSources,
      ArchGDAL, 
      NCDatasets

using SimpleSDMDataSources: Solar, Rainfall, VapourPressure, Temperature, 
      H09, H15, MinAve, MaxAve, WEATHER_LAYERS, AWAP_LAYERS

using GeoData: StandardIndices, readwindowed


const SSL = SimpleSDMDataSources


export CachedStack 

export InterpStack, InterpArray

export MinMaxInterpStack, MinMaxSpec

export Cosine, HypTan

include("cachedstack.jl")
include("interpstack.jl")
include("minmaxstack.jl")
include("load.jl")

end # module

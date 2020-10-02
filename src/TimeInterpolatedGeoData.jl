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
      H09, H15, MinAve, MaxAve, rasterlayers

using GeoData: StandardIndices, readwindowed


const SSL = SimpleSDMDataSources


export CachedStack 

export InterpStack, InterpArray

export MinMaxInterpStack, MinMaxSpec

export Cosine, HypTan

export geoarray, stack, series, interpseries, minmaxseries 

include("cachedstack.jl")
include("interpstack.jl")
include("minmaxstack.jl")
include("interpolations.jl")
include("load.jl")

end # module

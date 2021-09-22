using TimeInterpolatedGeoData, GeoData, Interpolations, Dates, Test

using TimeInterpolatedGeoData: calcfrac

include(joinpath(dirname(pathof(TimeInterpolatedGeoData)), "../test/test_utils.jl"))

@testset "InterpStack" begin

    @test calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1)) == 0
    @test calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 2)) == 1
    @test calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 12)) == 0.5
    @test calcfrac(1, 2, 1) == 0.0
    @test calcfrac(1, 2, 1.5) == 0.5
    @test calcfrac(1, 2, 2.0) == 1.0

    @test Interpolations.value_weights(Cosine(), 0)[1] == 1.0
    @test Interpolations.value_weights(Cosine(), 0.5)[1] ≈ 0.5
    @test Interpolations.value_weights(Cosine(), 1)[1] == 0.0
    @test Interpolations.value_weights(Cosine(), 0)[2] == 0.0
    @test Interpolations.value_weights(Cosine(), 0.5)[2] ≈ 0.5
    @test Interpolations.value_weights(Cosine(), 1)[2] == 1.0

    @test Interpolations.value_weights(HypTan(), 0)[1] == 1.0
    @test Interpolations.value_weights(HypTan(), 0.5)[1] ≈ 0.5
    @test Interpolations.value_weights(HypTan(), 1)[1] == 0.0
    @test Interpolations.value_weights(HypTan(), 0)[2] == 0.0
    @test Interpolations.value_weights(HypTan(), 0.5)[2] ≈ 0.5
    @test Interpolations.value_weights(HypTan(), 1)[2] == 1.0


    int = interpolate([1, 2, 3, 4], BSpline(Cosine()))
    @test int[2.0] == 2.0 
    @test int[2.25] < 2.75
    @test int[2.25] == -cos(π * 2.25) / 2 + 2 + 0.5 
    @test int[2.5] == 2.5
    @test int[2.75] > 2.75
    @test int[2.75] == -cos(π * 2.75) / 2 + 2 + 0.5 
    @test int[3.0] == 3.0

    # ncexamples = "https://www.unidata.ucar.edu/software/netcdf/examples/"
    # ncmulti = geturl(joinpath(ncexamples, "test_echam_spectral.nc"))

    # series = GeoSeries([NCDstack(ncmulti), NCDstack(ncmulti)], (Ti(DateTime(2020, 1):Month(1):DateTime(2020, 2)),); childtype=nothing)

    # iseries = interpseries(series, DateTime(2020, 1):Day(1):DateTime(2020, 2), (ameltfrac=BSpline(Cosine()),))
    # iseries2 = interpseries(series, DateTime(2020, 1):Day(1):DateTime(2020, 2), (ameltfrac=BSpline(Linear()),))
    # iseries[1][:ameltfrac].frac
    # iseries2[1][:ameltfrac]
    # iseries[2].stacks[1].loaded
    # iseries[10].stacks[2].loaded
    # iseries[10][:srafs, 1:10, 1:10, 1]
    # clear!(iseries[10])

    gdalpath = geturl("https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif")
    stack1 = GeoStack((a=gdalpath, b=gdalpath))
    stack2 = GeoStack((a=0x02 .* GeoArray(gdalpath), b=0x02 .* GeoArray(gdalpath)))
    istack = interpstack((stack1, stack2), (a=BSpline(Linear()), b=BSpline(Cosine())), 0.25)
    A = istack[:a]
    B = istack[:b]
    A = istack[:a, 1:10, 1:10, 1]
    B = istack[:b, 1:10, 1:10, 1]
    b = read(istack[:a]) 
    # using Plots
    # istack[:b] |> plot

    # test the values somehow here
end

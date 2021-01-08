using TimeInterpolatedGeoData, Test, Dates, Interpolations, GeoData, RasterDataSources

include(joinpath(dirname(pathof(TimeInterpolatedGeoData)), "../test/test_utils.jl"))

using TimeInterpolatedGeoData: chooseinterparrays, MinMaxFracs, calcfrac, Cosine

spec = MinMaxSpec((tmin=Hour(5), tmax=Hour(14)), BSpline(Cosine()))

@testset "chooseinterparrays choose the right arrays and calculates correct fractions" begin
    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 5))
    @test chooseinterparrays(spec, Day(1), tfrac) == ((tmin=1, tmax=1), 0.0)
    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 14))
    @test chooseinterparrays(spec, Day(1), tfrac) == ((tmax=1, tmin=2), 0.0)

    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 13, 59, 59, 99))
    @test chooseinterparrays(spec, Day(1), tfrac)[1] == (tmin=1, tmax=1)
    @test chooseinterparrays(spec, Day(1), tfrac)[2] ≈ 1.0 atol=1e4
    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 4, 59, 59, 99))
    @test chooseinterparrays(spec, Day(1), tfrac)[1] == (tmax=0, tmin=1)
    @test chooseinterparrays(spec, Day(1), tfrac)[2] ≈ 1.0 atol=1e4
    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 2, 4, 59, 59, 99))
    @test chooseinterparrays(spec, Day(1), tfrac)[1] == (tmax=1, tmin=2)
    @test chooseinterparrays(spec, Day(1), tfrac)[2] ≈ 1.0 atol=1e4

    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 8))
    @test chooseinterparrays(spec, Day(1), tfrac)[1] == (tmin=1, tmax=1) 
    @test chooseinterparrays(spec, Day(1), tfrac)[2] ≈ 1/3 
    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 11))
    @test chooseinterparrays(spec, Day(1), tfrac)[1] == (tmin=1, tmax=1) 
    @test chooseinterparrays(spec, Day(1), tfrac)[2] ≈ 2/3 

    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 19))
    @test chooseinterparrays(spec, Day(1), tfrac)[1] == (tmax=1, tmin=2) 
    @test chooseinterparrays(spec, Day(1), tfrac)[2] ≈ 1/3 
    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 24))
    @test chooseinterparrays(spec, Day(1), tfrac)[1] == (tmax=1, tmin=2) 
    @test chooseinterparrays(spec, Day(1), tfrac)[2] ≈ 2/3 
end

# In memory min-max
tempspec = MinMaxSpec((tmin=Hour(5), tmax=Hour(14)), BSpline(Linear()))

@testset "MinMaxFracs have correct fractions for correct min and max arrays" begin
    tempfracs = MinMaxFracs(tempspec, Day(1), 0/24)
    @test tempfracs.indices == (tmax=0, tmin=1)
    @test tempfracs.frac == 10/15
    @test tempfracs.interpmode == BSpline(Linear())

    tempfracs = MinMaxFracs(tempspec, Day(1), 5/24)
    @test tempfracs.indices == (tmin=1, tmax=1)
    @test tempfracs.frac == 0.0

    tempfracs = MinMaxFracs(tempspec, Day(1), 14/24)
    @test tempfracs.indices == (tmax=1, tmin=2)
    @test tempfracs.frac == 0.0

    # Mid point between 5 and 14 (9.30am)
    tempfracs = MinMaxFracs(tempspec, Day(1), 9.5/24)
    @test tempfracs.indices == (tmin=1, tmax=1)
    @test tempfracs.frac ≈ 0.5

    # Mid point between 14 and 29 (9.30pm)
    tempfracs = MinMaxFracs(tempspec, Day(1), 21.5/24)
    @test tempfracs.indices == (tmax=1, tmin=2)
    @test tempfracs.frac ≈ 0.5
end


@testset "MinMax" begin
    specs = (temp=tempspec,)
    min0 = [1.0 2.0; 3.0 4.0]
    dimz = (X, Y)
    dates = DateTime(2001):Day(1):DateTime(2001, 12)
    # tmin increases by 0.1 each day, tmax is j0.0 degreas above tmin for all dates
    stacks = [GeoStack((tmin=GeoArray(min0 .+ 0.1i, dimz), tmax=GeoArray(min0  .+ 10 .+ 0.1i, dimz))) for i in eachindex(dates)]
    ser = GeoSeries(stacks, Ti(dates))
    mmser = minmaxseries(ser, DateTime(2001, 2):Hour(1):DateTime(2001, 4), (), specs)

    mmstack = mmser[DateTime(2001, 2, 1)]
    @test mmstack.frac == 0.0
    mmstack = mmser[DateTime(2001, 2, 1, 5)]
    @test mmstack.frac == 5/24
    @test mmstack.specs[:temp].frac == 0.0
    @test mmstack.specs[:temp].indices == (tmin=1, tmax=1)
    @test mmstack[:temp] == ser[DateTime(2001, 2, 1)][:tmin]

    # Minimum at 5am
    @test mmser[DateTime(2001, 2, 1, 5)][:temp] == ser[DateTime(2001, 2, 1)][:tmin]
    # Maximum at 2pm
    @test mmser[DateTime(2001, 2, 1, 14)][:temp] == ser[DateTime(2001, 2, 1)][:tmax]
    # 1/3 of the way between max and min the next day
    @test mmser[DateTime(2001, 2, 1, 19)][:temp] ≈ ser[DateTime(2001, 2, 1)][:tmax] .* (2/3) + ser[DateTime(2001, 2, 2)][:tmin] .* (1/3)
end


if Sys.islinux()
    @testset "AWAP MinMax" begin
        layers = (:tmin, :tmax)
        dates = (DateTime(2018, 12, 31), DateTime(2020, 1, 1))
        getraster(AWAP, layers; date=dates)

        # Weather time-series
        ser = series(AWAP, layers; date=dates, window=(Band(1),))

        index(ser, Ti)
        ser[DateTime(2019, 1)][:tmin]
        ser[DateTime(2019, 1)][:tmax]

        mmseries = minmaxseries(ser, DateTime(2019, 1):Hour(1):DateTime(2019, 12, 31), (), specs)
        iseries = interpseries(ser, DateTime(2019, 1):Hour(1):DateTime(2019, 12, 31), (tmax=BSpline(Linear()), tmin=BSpline(Linear()),))
        awapstack = stack(AWAP; date=DateTime(2019, 1, 5)) 
        mmseries[DateTime(2019, 1, 5, 6)][:temp]
        mmseries[DateTime(2019, 1, 5, 1)][:temp]
        mmseries[DateTime(2019, 1, 5, 2)][:temp]
        mmseries[DateTime(2019, 1, 5, 3)][:temp]
    end
end

@testset "WorldClim Climate meandayminmaxseries" begin
    layers = (:tmin, :tmax)
    months = 1:12
    download_raster(WorldClim{Climate}, layers; month=months)

    ser = series(WorldClim{Climate}; layers=layers)
    ser = set(ser, Ti=(DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1)))
    # We use a DimensionalData dim instead of a vector of dates because a 
    # Dim can have a step size with irregular spaced data - here 1 Hour.
    dates = Ti(vec([d + h for h in Hour.(0:23), d in index(ser, Ti)]); 
        mode=Sampled(Ordered(), Regular(Hour(1)), Intervals())
    )
    tempspec = MinMaxSpec((tmin=Hour(5), tmax=Hour(14)), BSpline(Linear()))
    specs = (temp=tempspec,)

    mmseries = meandayminmaxseries(ser, dates, (), specs)
    mmseries[DateTime(2001, 1, 1, 6)][:temp]
    mmseries[DateTime(2001, 4, 1, 1)][:temp]
    mmseries[DateTime(2001, 9, 1, 13)][:temp]
    mmseries[DateTime(2001, 11, 1, 2)][:temp]
end

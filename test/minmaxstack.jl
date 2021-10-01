using TimeInterpolatedGeoData, Test, Dates, Interpolations, GeoData, RasterDataSources

include(joinpath(dirname(pathof(TimeInterpolatedGeoData)), "../test/test_utils.jl"))

using TimeInterpolatedGeoData: chooseinterparrays, MinMaxFrac, calcfrac, Cosine

interpolator = MinMaxInterpolator((tmin=Hour(5), tmax=Hour(14)), BSpline(Cosine()))

@testset "chooseinterparrays choose the right arrays and calculates correct fractions" begin
    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 5))
    @test chooseinterparrays(interpolator, Day(1), tfrac) == ((tmin=1, tmax=1), 0.0)
    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 14))
    @test chooseinterparrays(interpolator, Day(1), tfrac) == ((tmax=1, tmin=2), 0.0)

    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 13, 59, 59, 99))
    @test chooseinterparrays(interpolator, Day(1), tfrac)[1] == (tmin=1, tmax=1)
    @test chooseinterparrays(interpolator, Day(1), tfrac)[2] ≈ 1.0 atol=1e4
    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 4, 59, 59, 99))
    @test chooseinterparrays(interpolator, Day(1), tfrac)[1] == (tmax=0, tmin=1)
    @test chooseinterparrays(interpolator, Day(1), tfrac)[2] ≈ 1.0 atol=1e4
    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 2, 4, 59, 59, 99))
    @test chooseinterparrays(interpolator, Day(1), tfrac)[1] == (tmax=1, tmin=2)
    @test chooseinterparrays(interpolator, Day(1), tfrac)[2] ≈ 1.0 atol=1e4

    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 8))
    @test chooseinterparrays(interpolator, Day(1), tfrac)[1] == (tmin=1, tmax=1)
    @test chooseinterparrays(interpolator, Day(1), tfrac)[2] ≈ 1/3
    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 11))
    @test chooseinterparrays(interpolator, Day(1), tfrac)[1] == (tmin=1, tmax=1)
    @test chooseinterparrays(interpolator, Day(1), tfrac)[2] ≈ 2/3

    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 19))
    @test chooseinterparrays(interpolator, Day(1), tfrac)[1] == (tmax=1, tmin=2)
    @test chooseinterparrays(interpolator, Day(1), tfrac)[2] ≈ 1/3
    tfrac = calcfrac(DateTime(2001, 1, 1), DateTime(2001, 1, 2), DateTime(2001, 1, 1, 24))
    @test chooseinterparrays(interpolator, Day(1), tfrac)[1] == (tmax=1, tmin=2)
    @test chooseinterparrays(interpolator, Day(1), tfrac)[2] ≈ 2/3
end

# In memory min-max
tempinterpolator = MinMaxInterpolator((tmin=Hour(5), tmax=Hour(14)), BSpline(Linear()))

@testset "MinMaxFracs have correct fractions for correct min and max arrays" begin
    tempfracs = MinMaxFrac(tempinterpolator, Day(1), 0/24)
    @test tempfracs.indices == (tmax=0, tmin=1)
    @test tempfracs.frac == 10/15
    @test tempfracs.interptype == BSpline(Linear())

    tempfracs = MinMaxFrac(tempinterpolator, Day(1), 5/24)
    @test tempfracs.indices == (tmin=1, tmax=1)
    @test tempfracs.frac == 0.0

    tempfracs = MinMaxFrac(tempinterpolator, Day(1), 14/24)
    @test tempfracs.indices == (tmax=1, tmin=2)
    @test tempfracs.frac == 0.0

    # Mid point between 5 and 14 (9.30am)
    tempfracs = MinMaxFrac(tempinterpolator, Day(1), 9.5/24)
    @test tempfracs.indices == (tmin=1, tmax=1)
    @test tempfracs.frac ≈ 0.5

    # Mid point between 14 and 29 (9.30pm)
    tempfracs = MinMaxFrac(tempinterpolator, Day(1), 21.5/24)
    @test tempfracs.indices == (tmax=1, tmin=2)
    @test tempfracs.frac ≈ 0.5
end


@testset "MinMax" begin
    mm_interpolators = (temp=tempinterpolator,)
    min0 = [1.0 2.0; 3.0 4.0]
    dimz = (X, Y)
    dates = DateTime(2001):Day(1):DateTime(2001, 12)
    # tmin increases by 0.1 each day, tmax is j0.0 degreas above tmin for all dates
    stacks = [GeoStack((tmin=GeoArray(min0 .+ 0.1i, dimz), tmax=GeoArray(min0  .+ 10 .+ 0.1i, dimz))) for i in eachindex(dates)];
    ser = GeoSeries(stacks, Ti(dates))
    mmser = minmaxseries(ser, DateTime(2001, 2):Hour(1):DateTime(2001, 4); mm_interpolators)

    mmstack = mmser[At(DateTime(2001, 2, 1, 2))]
    @test mmstack[:temp].data.frac ≈ 1 + 12 / 15
    mmstack = mmser[At(DateTime(2001, 2, 1, 5))]
    @test mmstack[:temp].data.frac == 1.0
    @test mmstack[:temp] == ser[At(DateTime(2001, 2, 1))][:tmin]
    mmstack = mmser[At(DateTime(2001, 2, 1, 14))]
    @test mmstack[:temp].data.frac == 1.0
    @test mmstack[:temp] == ser[At(DateTime(2001, 2, 1))][:tmax]

    # Minimum at 5am
    @test mmser[At(DateTime(2001, 2, 1, 5))][:temp] == ser[At(DateTime(2001, 2, 1))][:tmin]
    # Maximum at 2pm
    @test mmser[At(DateTime(2001, 2, 1, 14))][:temp] == ser[At(DateTime(2001, 2, 1))][:tmax]
    # 1/3 of the way between max and min the next day
    @test mmser[At(DateTime(2001, 2, 1, 19))][:temp] ≈
        ser[At(DateTime(2001, 2, 1))][:tmax] .* (2/3) + ser[At(DateTime(2001, 2, 2))][:tmin] .* (1/3)
end


if Sys.islinux()
    @testset "AWAP MinMax" begin
        layers = (:tmin, :tmax, :rainfall)
        serdates = (DateTime(2019, 12, 31), DateTime(2020, 2, 1))
        dates = DateTime(2020, 1, 1):Hour(1):DateTime(2020, 1, 31, 23)

        # Weather time-series
        ser = GeoSeries(AWAP, layers; date=serdates, resize=nothing)

        mm_interpolators = (;
            temp=MinMaxInterpolator((tmin=Hour(5), tmax=Hour(14)), BSpline(Linear()))
        )
        mmseries = minmaxseries(ser, dates; mm_interpolators=mm_interpolators)
        @test length(mmseries) == 744
        @test keys(mmseries[1]) === (:temp, :rainfall)
        iseries = interpseries(ser, dates; 
             interptypes=(tmax=BSpline(Linear()), tmin=BSpline(Linear()),)
        )
        @test length(iseries) == 744
        @test keys(iseries[1]) === (:tmin, :tmax, :rainfall)

        @test all(iseries[At(DateTime(2020, 1, 1, 0))][:tmin] .==
            mmseries[At(DateTime(2020, 1, 1, 5))][:temp] .== 
            ser[At(DateTime(2020, 1))][:tmin]
        )
        @test all(iseries[At(DateTime(2020, 1, 1, 0))][:tmax] .==
            mmseries[At(DateTime(2020, 1, 1, 14))][:temp] .== 
            ser[At(DateTime(2020, 1))][:tmax]
        )
        @test all(iseries[At(DateTime(2020, 1, 1, 0))][:tmax] .==
            mmseries[At(DateTime(2020, 1, 1, 14))][:temp] .== 
            ser[At(DateTime(2020, 1))][:tmax]
        )
        @test all(mmseries[At(DateTime(2020, 1, 1, 7))][:temp] .≈ 
            ser[At(DateTime(2020, 1))][:tmin] .* 7 / 9 .+ 
            ser[At(DateTime(2020, 1))][:tmax] .* 2 / 9
        )
    end
end

@testset "WorldClim Climate meanday_minmaxseries" begin
    layers = (:tmin, :tmax)
    months = 1:12
    getraster(WorldClim{Climate}, layers; month=months)

    ser = series(WorldClim{Climate}, layers; month=Jan:Dec)
    ser = set(ser, :month => Ti(DateTime(2001, 1, 1):Month(1):DateTime(2001, 12, 1)))
    # We use a DimensionalData dim instead of a vector of dates because a
    # Dim can have a step size with irregular spaced data - here 1 Hour.
    dates = [d + h for d in index(ser, Ti) for h in Hour.(0:23)]
    
    tempinterpolator = MinMaxInterpolator((tmin=Hour(5), tmax=Hour(14)), BSpline(Linear()))
    mm_interpolators = (temp=tempinterpolator,)
    mmseries = meanday_minmaxseries(ser, dates; step=Hour(1), mm_interpolators)
    @test keys(mmseries[1]) == (:temp,)
    @test mmseries[At(DateTime(2001, 1, 1, 5))][:temp] == ser[At(DateTime(2001, 1, 1))][:tmin]
    @test mmseries[At(DateTime(2001, 4, 1, 14))][:temp] == ser[At(DateTime(2001, 4, 1))][:tmax]
    
    # Animation
    #
    # using Plots
    # anim = @animate for m in 1:12, h in 1:23 
    #     plot(mmseries[At(DateTime(2001, m, 1, h))][:temp]; title=string("2001-", m, "-1-", h), clim=(0, 50)) |> display
    # end
    # gif(anim, "minmax.gif", fps = 12)


    # Compare mean of :temp to :tavg
    #
    # using Statistics
    # monthblocks = map(1:12) do m
    #     monthser = read(mmseries[Ti(Where(t -> month(t) == m))])
    #     catted = cat(monthser...; dims=dims(monthser, Ti))
    #     mean(catted; dims=Ti)[Band(1), Ti(1)]
    # end
    # tavg_ser = series(WorldClim{Climate}, :tavg; month=Jan:Dec)
    # tavg_ser = map(x -> view(x, Band(1)), tavg_ser)
    
    # using Plots
    # masklayer = missingmask(tavg_ser[1])
    # anim = @animate for m in 1:12
    #     title = "Month: " * monthname(m)
    #     A = monthblocks[m][:temp] .- tavg_ser[m] .* masklayer
    #     plot(A; title) 
    # end
    # gif(anim, "diff_minmax_tavg.gif", fps = 4)
    # spoiler: its surprisingly similar
end

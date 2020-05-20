using Pkg.Artifacts
using Test
using MWAjl

testdata = artifact"testdata"

@testset "SkyModelParser" begin
    open(joinpath(testdata, "data/model.txt")) do f
        @test length(parse_model(f)) == 200
    end

    @test MWAjl.hms2rad("1h3m14.5s") ≈ 0.27594382694552017
    @test MWAjl.hms2rad("00h05m59.1299s") ≈ 0.026116663322324932
    @test MWAjl.hms2rad("-0h05m59.1299s") ≈ -0.026116663322324932
    @test MWAjl.dms2rad("-14d37m20s") ≈ -0.25520592173605977
end

@testset "SkyModel" begin
    @testset "Measurement flux interpolation" begin
        ms = MWAjl.Measurements(MWAjl.Measurement[
            MWAjl.Measurement(100E6, Float64[1, 0, 0, 0]),
            MWAjl.Measurement(300E6, Float64[0.1, 0, 0, 0]),
        ])

        @test MWAjl.stokes(ms, 100E6) ≈ Float64[1, 0, 0, 0]
        @test MWAjl.stokes(ms, 200E6) ≈ Float64[0.23392155723884647, 0, 0, 0]
        @test MWAjl.stokes(ms, 800E6) ≈ Float64[0.012800022683622247, 0, 0, 0]

        # Edge case where one measurement is zero
        ms = MWAjl.Measurements(MWAjl.Measurement[
            MWAjl.Measurement(100E6, Float64[1, 0, 0, 0]),
            MWAjl.Measurement(300E6, Float64[0, 0, 0, 0]),
        ])
        @test MWAjl.stokes(ms, 300E6) ≈ Float64[0, 0, 0, 0]
        @test MWAjl.stokes(ms, 400E6) ≈ Float64[1, 0, 0, 0]

        # Edge case: just one measurement
        ms = MWAjl.Measurements(MWAjl.Measurement[
            MWAjl.Measurement(100E6, Float64[1, 0, 0, 0]),
        ])
        @test MWAjl.stokes(ms, 300E6) ≈ Float64[1, 0, 0, 0]

        # Edge case: just one measurement, and it's negative!
        ms = MWAjl.Measurements(MWAjl.Measurement[
            MWAjl.Measurement(100E6, Float64[-1, 0, 0, 0]),
        ])
        @test MWAjl.stokes(ms, 100E6) ≈ Float64[-1, 0, 0, 0]
        @test MWAjl.stokes(ms, 200E6) ≈ Float64[-1, 0, 0, 0]

    end

    @testset "SED flux interpolation" begin
        sed = MWAjl.SED(154E6, [1, 0, 0, 0], [-1])
        @test MWAjl.stokes(sed, 180E6) ≈ [0.8555555555555556, 0, 0, 0]

        sed = MWAjl.SED(154E6, [1, 0, 0, 0], [-1, 0.001, 0.00003])
        @test MWAjl.stokes(sed, 180E6) ≈ [0.855576475194691, 0, 0, 0]
    end
end
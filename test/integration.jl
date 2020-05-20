using Pkg.Artifacts
using Test

include(joinpath(@__DIR__, "../bin/calibrate.jl"))
const testdata = artifact"testdata"

args = Dict(
    "model" => joinpath(testdata, "data/model.txt"),
    "nbatch" => 200,
    "tolerance" => [1E-5, 1E-6],
    "mset" => joinpath(testdata, "data/1217094992-ch8-19.ms"),
    "solution" => "/dev/null",
    "chanwidth" => 1,
    "timewidth" => 0,
    "minuv" => 0.0,
    "maxuv" => 9E99,
    "max-iterations" => 50,
    "apply-beam" => false,
    "datacolumn" => "DATA",
    "modelcolumn" => "MODEL_DATA",
    "debug" => false,
    "verbose" => true,
    "gpu" => false,
)

@testset "Integration tests" begin
    @testset "CPU Prediction + Calibration" begin
        @test main(args)
    end

    @testset "CPU Prediction + Calibration + Beam" begin
        let args = merge(args, Dict("apply-beam" => true))
            @test main(args)
        end
    end

    @testset "GPU Prediction + Calibration" begin
        let args = merge(args, Dict("gpu" => true))
            @test main(args)
        end
    end

    @testset "GPU Prediction + Calibration + Beam" begin
        let args = merge(args, Dict("gpu" => true, "apply-beam" => true))
            @test main(args)
        end
    end

    @testset "With timewith = 10" begin
        let args = merge(args, Dict("gpu" => true, "timewidth" => 10))
            @test main(args)
        end
    end

    @testset "With chanwith = 3" begin
        let args = merge(args, Dict("gpu" => true, "chanwidth" => 3))
            @test main(args)
        end
    end

    @testset "With minuv = 30" begin
        let args = merge(args, Dict("gpu" => true, "minuv" => 30.))
            @test main(args)
        end
    end
end

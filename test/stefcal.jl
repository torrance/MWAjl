using MWAjl
using NPZ
using Pkg.Artifacts
using Statistics: mean
using Test

testdata = artifact"testdata"

@testset "Simple Calibration" begin
    nants = 128
    ntimesteps = 30
    nchans = 4
    nrows = Int((nants * (nants - 1)) / 2 * ntimesteps)

    for tid in 1:2
        data = zeros(ComplexF32, 2, 2, nchans, nrows)
        model = zeros(ComplexF32, 2, 2, nchans, nrows)
        ants1 = zeros(Int, nrows)
        ants2 = zeros(Int, nrows)
        weights = randn(Float32, 4, nchans, nrows)
        weights .= abs.(1 .+ weights)

        # Initialize data
        data[1, 1, :, :] .= 1
        data[2, 2, :, :] .= 1
        model[1, 1, :, :] .= 1
        model[2, 2, :, :] .= 1

        # if tid == 2
        #     # Fuck up data
        #     data[:, :, :, 1:2:end] .*= 100
        #     weights[:, :, :, 1:2:end] .= 0.005
        # end

        # Create true jones
        jones = (reshape([1 0; 0 1], 2, 2, 1) .+ 0.2 * randn(Float32, 2, 2, 128)) .* exp.(2im * π * rand(Float32, 2, 2, 128))

        # Uncalibrate data
        row = 1
        for timestep in 1:ntimesteps, ant1 in 1:nants, ant2 in 1:nants
            if ant1 >= ant2
                continue
            end
            ants1[row] = ant1
            ants2[row] = ant2
            for chan in 1:nchans
                data[:, :, chan, row] = jones[:, :, ant1] * data[:, :, chan, row] * jones[:, :, ant2]'
            end
            row += 1
        end

        # Find solution
        jones1 = similar(jones)
        jones1 .= 0
        jones1[1, 1, :] .= 1
        jones1[2, 2, :] .= 1
        jones1 = Matrix2x2toArray4(jones1)
        calibrate!(jones1, Matrix2x2toArray4(data), Matrix2x2toArray4(model), weights, ants1, ants2, 50, 1E-5, 1E-8)

        # Correct data with new Jones solution
        jones1 = Array4toMatrix2x2(jones1)
        for antid in axes(jones1, 3)
            jones1[:, :, antid] = inv(jones1[:, :, antid])
        end
        corrected = similar(data)
        for row in axes(data, 4), chan in axes(data, 3)
            corrected[:, :, chan, row] = jones1[:, :, ants1[row]] * data[:, :, chan, row] * jones1[:, :, ants2[row]]'
        end

        @test isapprox(model, corrected, nans=true, atol=1E-4)
    end
end

@testset "Real data with partial skymodel" begin
    ants1 = npzread(joinpath(testdata, "data/1217094992/1217094992-ANTENNA1.npy")) .+ 1
    ants2 = npzread(joinpath(testdata, "data/1217094992/1217094992-ANTENNA2.npy")) .+ 1
    auto = ants1 .== ants2

    ants1 = ants1[.~auto]
    ants2 = ants2[.~auto]
    data = permutedims(npzread(joinpath(testdata, "data/1217094992/1217094992-DATA-ch2-3.npy"))[.~auto, :, :], [3, 2, 1])
    model = permutedims(npzread(joinpath(testdata, "data/1217094992/1217094992-MODEL_DATA-ch2-3.npy"))[.~auto, :, :], [3, 2, 1])
    flag = permutedims(npzread(joinpath(testdata, "data/1217094992/1217094992-FLAG-ch2-3.npy"))[.~auto, :, :], [3, 2, 1])
    weights = permutedims(npzread(joinpath(testdata, "data/1217094992/1217094992-WEIGHT_SPECTRUM-ch2-3.npy"))[.~auto, :, :], [3, 2, 1])

    sanitize!(data, model, flag)

    jones = zeros(ComplexF64, 4, 128)
    jones[1, :] .= 1
    jones[4, :] .= 1
    success, iterations = calibrate!(jones, data, model, weights, ants1, ants2, 100, 1E-5, 1E-8)
    @test success
end

# TODO: Test where all solutions fail

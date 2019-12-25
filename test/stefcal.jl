using MWAjl
using Statistics: mean
using Test

@testset "Simple Calibration" begin
    nrows = 140_000
    chans = 2

    ants1 = rand(1:128, nrows)
    ants2 = rand(1:128, nrows)

    for tid in 1:2
        data = zeros(ComplexF32, 2, 2, chans, nrows)
        model = zeros(ComplexF32, 2, 2, chans, nrows)
        weights = randn(Float32, size(data)...)
        weights .= abs.(1 .+ weights)

        # Initialize data
        data[1, 1, :, :] .= 1
        data[2, 2, :, :] .= 1
        model[1, 1, :, :] .= 1
        model[2, 2, :, :] .= 1

        if tid == 2
            # Fuck up data
            data[:, :, :, 1:2:end] .*= 100
            weights[:, :, :, 1:2:end] .= 0.005
        end

        # Create true jones
        jones = (reshape([1 0; 0 1], 2, 2, 1) .+ 0.2 * randn(Float32, 2, 2, 128)) .* exp.(2im * π * rand(Float32, 2, 2, 128))

        # Uncalibrate data
        for row in axes(data, 4), chan in axes(data, 3)
            data[:, :, chan, row] = jones[:, :, ants1[row]] * data[:, :, chan, row] * jones[:, :, ants2[row]]'
        end

        # Find solution
        jones1 = similar(jones)
        jones1 .= 0
        jones1[1, 1, :] .= 1
        jones1[2, 2, :] .= 1
        calibrate!(jones1, data, model, weights, ants1, ants2, 1E-5, 1E-8)

        corrected = similar(data)
        for row in axes(data, 4), chan in axes(data, 3)
            corrected[:, :, chan, row] = inv(jones1[:, :, ants1[row]]) * data[:, :, chan, row] * inv(jones1[:, :, ants2[row]]')
        end

        @test isapprox(model, corrected, nans=true, atol=1E-4)
    end
end

# TODO: Test where all solutions fail

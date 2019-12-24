using MWAjl
using Statistics: mean
using Test

@testset "Simple Calibration" begin
    nrows = 231_000
    chans = 3

    ants1 = rand(1:128, nrows)
    ants2 = rand(1:128, nrows)

    for _ in 1:5
        data = zeros(ComplexF64, 2, 2, chans, nrows)
        model = zeros(ComplexF64, 2, 2, chans, nrows)

        # Initialize data
        data[1, 1, :, :] .= 1
        data[2, 2, :, :] .= 1
        model[1, 1, :, :] .= 1
        model[2, 2, :, :] .= 1

        # Create true jones
        jones = (reshape([1 0; 0 1], 2, 2, 1) .+ 0.2 * randn(Float64, 2, 2, 128)) .* exp.(2im * π * rand(Float64, 2, 2, 128))

        # Uncalibrate data
        for row in axes(data, 4), chan in axes(data, 3)
            data[:, :, chan, row] = jones[:, :, ants1[row]] * data[:, :, chan, row] * jones[:, :, ants2[row]]'
        end

        jones1 = calibrate(data, model, ants1, ants2)

        corrected = similar(data)
        for row in axes(data, 4), chan in axes(data, 3)
            corrected[:, :, chan, row] = inv(jones1[:, :, ants1[row]]) * data[:, :, chan, row] * inv(jones1[:, :, ants2[row]]')
        end

        @test isapprox(model, corrected, nans=true, atol=1E-3)
    end
end

# TODO: Test where all solutions fail

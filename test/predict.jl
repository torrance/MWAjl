using Logging
using MWAjl
using Test

global_logger(Logging.ConsoleLogger(stderr, Logging.Debug))

@testset "Prediction" begin
    nrows = 500_000
    nchans = 4

    uvws = rand(3, nrows)
    times = sort(rand(1.0:10.0, nrows))
    freqs = collect(range(140E6, length=nchans, stop=160E6))
    lambdas = 299792458 ./ freqs
    pos0 = Position(0, 0)

    # Single source at phase center
    comps = Component[Component(
        Position(0, 0),
        "point",
        MWAjl.SED(154E6, [1, 0, 0, 0], [0]),
    )]
    model = predict(uvws, times, lambdas, comps, nothing, pos0)
    @test all(model[[true, false, false, true], :, :] .≈ 1)
    @test all(model[[false, true, true, false], :, :] .≈ 0)

    # Single source off phase center
    comps = Component[Component(
        Position(0.2, 0.2),
        "point",
        MWAjl.SED(154E6, [1, 0, 0, 0], [0]),
    )]
    model = predict(uvws, times, lambdas, comps, nothing, pos0)

    @test all(model[[false, true, true, false], :, :] .≈ 0)
    u, v, w = uvws[:, 1] ./ lambdas[1]
    l, m, n = lmn(comps[1], pos0)
    @test model[1, 1, 1] ≈ exp(2im * π * (u*l + v*m + w*(n - 1)))

    # Two sources off phase center
    comps = Component[
        Component(
            Position(0.2, 0.2),
            "point",
            MWAjl.SED(154E6, [1, 0, 0, 0], [0]),
        ),
        Component(
            Position(0.1, -0.15),
            "point",
            MWAjl.SED(154E6, [1, 0, 0, 0], [0]),
        ),
    ]
    model = predict(uvws, times, lambdas, comps, nothing, pos0)

    @test all(model[[false, true, true, false], :, :] .≈ 0)
    vis = 0
    u, v, w = uvws[:, 1] ./ lambdas[1]
    for i in 1:2
        l, m, n = lmn(comps[i], pos0)
        vis += exp(2im * π * (u*l + v*m + w*(n - 1)))
    end
    @test model[1, 1, 1] ≈ vis
end

@testset "GPU Prediction" begin
    nrows = 500_000
    nchans = 100

    uvws = rand(3, nrows)
    times = sort(rand(1.0:10.0, nrows))
    freqs = collect(range(140E6, length=nchans, stop=160E6))
    lambdas = 299792458 ./ freqs
    pos0 = Position(0, 0)

    comps = Component[]
    open(string(@__DIR__, "/data/model.txt")) do f
        sources = parse_model(f)
        for source in sources, comp in source.components
            push!(comps, comp)
        end
    end

    gpumodel = predict(uvws, times, lambdas, comps, nothing, pos0, gpu=true)
    cpumodel = predict(uvws, times, lambdas, comps, nothing, pos0)

    @test all(cpumodel .≈ gpumodel)
end
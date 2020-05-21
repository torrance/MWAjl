using CUDAdrv
using CUDAnative
using StaticArrays

struct Gaussian <: Component
    position::Position
    spectrum::Spectrum
    major::Float32
    minor::Float32
    pa::Float32
end

function predictgaussian!(
        model::Array{ComplexF32, 3},
        uvws::Array{Float32, 2},
        lambdas::Array{Float32, 1},
        timesteps::Array{Int, 1},
        lmns::Array{Float32, 2},
        fluxes::Array{ComplexF32, 4},
        majors::Array{Float32, 1},
        minors::Array{Float32, 1},
        pas::Array{Float32, 1},
    )

    # Threads.@threads doesn't work when already inside a thread (Julia v1.4)
    # so here we are doing it ourselves.
    # Distribute rows amongst available threads
    N = Threads.nthreads()
    rows = size(uvws, 2)
    step = ceil(Int, rows / N)

    tasks = Task[]
    elapsed = Base.@elapsed for idx1 in range(1, rows, step=step)
        idx2 = min(idx1 + step - 1, rows)

        task = Threads.@spawn _predictgaussian(view(model, :, :, idx1:idx2), view(uvws, :, idx1:idx2), $lambdas, view(timesteps, idx1:idx2), $lmns, $fluxes, $majors, $minors, $pas)
        push!(tasks, task)
    end

    # Wait for the tasks to finish
    for task in tasks
        wait(task)
    end
end

@fastmath @inline @inbounds function _predictgaussian(
    model::SubArray{ComplexF32, 3},
    uvws::SubArray{Float32, 2},
    lambdas::Array{Float32, 1},
    timesteps::SubArray{Int, 1},
    lmns::Array{Float32, 2},
    fluxes::Array{ComplexF32, 4},
    majors::Array{Float32, 1},
    minors::Array{Float32, 1},
    pas::Array{Float32, 1},
)
    cospas, sinpas = cos.(pas), sin.(pas)

    local timestep::Int
    local u::Float32, v::Float32, w::Float32
    local envelope::Float32
    local phase::ComplexF32

    for row in axes(model, 3)
        timestep = timesteps[row]
        for chan in axes(model, 2)
            u, v, w = uvws[:, row] / lambdas[chan]

            for compidx in axes(fluxes, 2)
                # Calculate Gaussian envelope, but short circuit to 1 for point sources
                if majors[compidx] == 0 || minors[compidx] == 0
                    envelope = 1.
                else
                    envelope = exp(-2 * π^2 * (
                        majors[compidx]^2 * (u * cospas[compidx] - v * sinpas[compidx])^2 +
                        minors[compidx]^2 * (u * sinpas[compidx] + v * cospas[compidx])^2
                    ))
                end
                phase = exp(
                    2im * π * (
                        u * lmns[1, compidx] +
                        v * lmns[2, compidx] +
                        w * (lmns[3, compidx] - 1)
                    )
                )
                for pol in 1:4
                    model[pol, chan, row] += fluxes[pol, compidx, chan, timestep] * envelope * phase
                end
            end
        end
    end
end

@inbounds @views @fastmath function predictgaussian!(
        model::CuDeviceArray{ComplexF32, 3},
        uvws::CuDeviceArray{Float32, 2},
        lambdas::CuDeviceArray{Float32, 1},
        timesteps::CuDeviceArray{Int, 1},
        lmns::CuDeviceArray{Float32, 2},
        fluxes::CuDeviceArray{ComplexF32, 4},
        majors::CuDeviceArray{Float32, 1},
        minors::CuDeviceArray{Float32, 1},
        pas::CuDeviceArray{Float32, 1},
    )

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x

    # Prefetch regularly used data into the (faster) shared memory.
    # This is limited in size, so only lmns and lambdas can reliably fit.
    # However, there's only a very marginal speed increase from doing
    # this, presumably because we are dominated by the cost of the complex
    # exponential function. Nonetheless, we leave it here as an example and
    # reminder that this was tried.
    lmns_shm = @cuDynamicSharedMem(Float32, size(lmns))
    for i in threadIdx().x:blockDim().x:length(lmns)
        lmns_shm[i] = lmns[i]
    end

    lambdas_shm = @cuDynamicSharedMem(Float32, size(lambdas), offset=sizeof(Float32) * length(lmns))
    for i in threadIdx().x:blockDim().x:length(lambdas)
        lambdas_shm[i] = lambdas[i]
    end

    # Wait for all threads to finish setting the shared memory cache, before we
    # read from it, so as to avoid race conditions.
    sync_threads()

    # We allocate a static vector locally per thread to do inplace
    # addition, so that we write out to model just once per index.
    # This has a significant performance improvement.
    cell = MVector{4, ComplexF32}(undef)
    timestep = 0::Int
    envelope = Float32(0)

    for row in index:stride:size(model, 3)
        timestep = timesteps[row]

        for chan in axes(model, 2)
            fill!(cell, 0)
            u, v, w = uvws[1, row] / lambdas_shm[chan], uvws[2, row] / lambdas_shm[chan], uvws[3, row] / lambdas_shm[chan]
            for compidx in axes(fluxes, 2)
                # Calculate Gaussian envelope, but short circuit to 1 for point sources
                if majors[compidx] == 0 || minors[compidx] == 0
                    envelope = 1
                else
                    cospa = CUDAnative.cos_fast(pas[compidx])
                    sinpa = CUDAnative.sin_fast(pas[compidx])
                    envelope = CUDAnative.exp(-2 * π^2 * (
                        majors[compidx]^2 * (u * cospa - v * sinpa)^2 +
                        minors[compidx]^2 * (u * sinpa + v * cospa)^2
                    ))
                end
                phase = CUDAnative.exp(
                    2im * π * (
                        u * lmns_shm[1, compidx] +
                        v * lmns_shm[2, compidx] +
                        w * (lmns_shm[3, compidx] - 1)
                    )
                )
                for pol in 1:4
                    cell[pol] += fluxes[pol, compidx, chan, timestep] * envelope * phase
                end
            end

            for pol in 1:4
                model[pol, chan, row] = cell[pol]
            end
        end
    end
end
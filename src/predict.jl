using CuArrays
using CUDAdrv
using CUDAnative
using StaticArrays

function predict(
        uvws::Array{T, 2},
        times::Array{T, 1},
        lambdas::Array{T, 1},
        comps::Array{Component, 1},
        beam::Union{Beam, Nothing},
        pos0::Position;
        gpu::Bool = false,
    ) where {T <: AbstractFloat}

    freqs = 299792458 ./ lambdas

    # Construct an time index into each row
    elapsed = Base.@elapsed begin
        unique_times = unique(times)
        timeidxs = convert(Array{Int}, indexin(times, unique_times))
    end
    @debug "Calculated time index elapsed $elapsed"

    # Precalculate l, m, n coordinates
    # TODO: This only needs to be done once
    elapsed = Base.@elapsed begin
        lmns = zeros(3, length(comps))
        for (compidx, comp) in enumerate(comps)
            lmns[:, compidx] .= lmn(comp, pos0)
        end
    end
    @debug "Calculated lmns elapsed $elapsed"

    # Precalculate flux matrices for each source, for each frequency
    elapsed = Base.@elapsed begin
        # Fluxes = [pol, comps, chans, times]
        fluxes = zeros(Complex{T}, 4, length(comps), length(lambdas), length(unique_times))
        for (chan, freq) in enumerate(freqs), (compidx, comp) in enumerate(comps)
            fluxes[:, compidx, chan, :] .= instrumental(comp, freq)[:, :]
        end
    end
    @debug "Calculated flux matrices for sources elapsed $elapsed"
    @debug "Absolute brightness at low frequency: $(sum(fluxes[1, :, 1, 1]))"
    @debug "Absolute brightness at high frequency: $(sum(fluxes[1, :, end, 1]))"

    # Apply beam correction
    if beam !== nothing
        # Calculate beam values for all timesteps (in one batch)
        # First, calculate apparent positions of sources
        # TODO: cache this data and only calculate once
        alts, azs = Float64[], Float64[]
        elapsed = Base.@elapsed for (timeidx, time) in enumerate(unique_times)
            # Calculate apparent sky coordinates for this timestep
            # Measurement set columns are in modified Julian date, but in *seconds*.
            # We hardcode the location of MWA for now.
            frame = Frame(time / (24 * 3600), 2.0362898433426406, -0.46606084551737853)
            for comp in comps
                alt, az = radec2altaz(comp.position, frame)
                push!(alts, alt)
                push!(azs, az)
            end
        end
        @debug "Calculated AltAz positions of $(length(comps)) at $(length(unique_times)) timesteps, elapsed $elapsed"

        # Get beam Jones matrix, and correct model fluxes to apparent fluxes
        tmp = zeros(Complex{T}, 4)
        for (chan, freq) in enumerate(freqs)
            elapsed = Base.@elapsed begin
                jones = reshape(
                    beamjones(beam, freq, alts, azs),
                    4, length(comps), length(unique_times)
                )
            end
            @debug "Retrieved beam Jones values, elapsed $elapsed"
            # apparent = J A J^H
            elapsed = Base.@elapsed @views for timeidx in axes(jones, 3), compidx in axes(jones, 2)
                AxBH!(tmp, fluxes[:, compidx, chan, timeidx], jones[:, compidx, timeidx])
                AxB!(fluxes[:, compidx, chan, timeidx], jones[:, compidx, timeidx], tmp)
            end
            @debug "Applied beam correction to fluxes, elapsed $elapsed"
        end
    end

    # Allocate model array
    if gpu
        model_d = CuArrays.fill(ComplexF32(0), 4, length(lambdas), length(times))
        uvws_d, lambdas_d, timeidxs_d, lmns_d, fluxes_d = CuArray(convert(Array{Float32}, uvws)), CuArray(convert(Array{Float32}, lambdas)), CuArray(timeidxs), CuArray(convert(Array{Float32}, lmns)), CuArray(convert(Array{ComplexF32}, fluxes))
        elapsed = Base.@elapsed CuArrays.@sync begin
            @cuda blocks=256 threads=256 shmem=sizeof(Float32)*length(lmns) + sizeof(Float32)*length(lambdas) predictionloop!(model_d, uvws_d, lambdas_d, timeidxs_d, lmns_d, fluxes_d)
        end
        @debug "Predicted model (GPU) elapsed $elapsed"

        return Array(model_d)
    else
        model = zeros(ComplexF32, 4, length(lambdas), length(times))
        elapsed = Base.@elapsed predictionloop!(model, uvws, lambdas, timeidxs, lmns, fluxes)
        @debug "Predicted model (CPU) elapsed $elapsed"

        return model
    end
end

@inbounds @views function predictionloop!(model::Array{ComplexF32, 3}, uvws, lambdas, timeidxs, lmns, fluxes)
    Threads.@threads for row in axes(model, 3); for chan in axes(model, 2)
        for compidx in axes(fluxes, 2)
            phase = exp(
                2im * π * (
                    uvws[1, row] * lmns[1, compidx] +
                    uvws[2, row] * lmns[2, compidx] +
                    uvws[3, row] * (lmns[3, compidx] - 1)
                ) / lambdas[chan]
            )
            for pol in 1:4
                model[pol, chan, row] += fluxes[pol, compidx, chan, timeidxs[row]] * phase
            end
        end
    end; end
end

@inbounds @views @fastmath function predictionloop!(
        model::CuDeviceArray{ComplexF32, 3},
        uvws::CuDeviceArray{Float32, 2},
        lambdas::CuDeviceArray{Float32, 1},
        timeidxs::CuDeviceArray{Int, 1},
        lmns::CuDeviceArray{Float32, 2},
        fluxes::CuDeviceArray{ComplexF32, 4},
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

    for row in index:stride:size(model, 3)
        for chan in axes(model, 2)
            fill!(cell, 0)
            for compidx in axes(fluxes, 2)
                phase = CUDAnative.exp(
                    2im * π * (
                        uvws[1, row] * lmns_shm[1, compidx] +
                        uvws[2, row] * lmns_shm[2, compidx] +
                        uvws[3, row] * (lmns_shm[3, compidx] - 1)
                    ) / lambdas_shm[chan]
                )
                for pol in 1:4
                    cell[pol] += fluxes[pol, compidx, chan, timeidxs[row]] * phase
                end
            end

            for pol in 1:4
                model[pol, chan, row] = cell[pol]
            end
        end
    end
end
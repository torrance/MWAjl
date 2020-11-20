using CuArrays

function predict(
        uvws::Array{T, 2},
        times::Array{T, 1},
        freqs::Array{T, 1},
        comps::Array{Component, 1},
        beam::Union{Beam, Nothing},
        pos0::Position;
        gpu::Bool = false,
    )::Union{Array{ComplexF32, 3}, CuArray{ComplexF32, 3}} where {T <: AbstractFloat}

    lambdas = 299792458 ./ freqs

    # Construct an time index into each row
    elapsed = Base.@elapsed begin
        unique_times = unique(times)
        timesteps = convert(Array{Int}, indexin(times, unique_times))
    end
    @debug "Calculated time index elapsed $elapsed"

    # Precalculate l, m, n coordinates
    # TODO: This only needs to be done once
    elapsed = Base.@elapsed begin
        lmns = zeros(Float32, 3, length(comps))
        for (compidx, comp) in enumerate(comps)
            lmns[:, compidx] .= lmn(comp, pos0)
        end
    end
    @debug "Calculated lmns elapsed $elapsed"

    # Precalculate flux matrices for each source, for each frequency
    elapsed = Base.@elapsed begin
        # Fluxes = [pol, comps, chans, times]
        fluxes = zeros(ComplexF32, 4, length(comps), length(lambdas), length(unique_times))
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
        elapsed = Base.@elapsed for time in unique_times
            # Calculate apparent sky coordinates for this timestep
            # Measurement set columns are in modified Julian date, but in *seconds*.
            # We hardcode the location of MWA for now.
            frame = Frame(Float64(time / (24 * 3600)), 2.0362898433426406, -0.46606084551737853)
            for comp in comps
                alt, az = radec2altaz(comp.position, frame)
                push!(alts, alt)
                push!(azs, az)
            end
        end
        @debug "Calculated AltAz positions of $(length(comps)) at $(length(unique_times)) timesteps, elapsed $elapsed"

        # Get beam Jones matrix, and correct model fluxes to apparent fluxes
        let tmp = zeros(ComplexF32, 4)
            elapsed = Base.@elapsed for (chan, freq) in enumerate(freqs)
                jones = reshape(
                    calc_jones(beam, azs, pi/2 .- alts, UInt32(freq)),
                    4, length(comps), length(unique_times)
                )

                # apparent = J A J^H
                @views for timeidx in axes(jones, 3), compidx in axes(jones, 2)
                    AxBH!(tmp, fluxes[:, compidx, chan, timeidx], jones[:, compidx, timeidx])
                    AxB!(fluxes[:, compidx, chan, timeidx], jones[:, compidx, timeidx], tmp)
                end
            end
            @debug "Calculated and applied beam corrections, elapsed $elapsed"
        end
    end

    # Extract Gaussian shape parameters
    # Convert FWHM major and minors to sigmas
    majors = Float32[comp.major / (2 * sqrt(2 * log(2))) for comp in comps]
    minors = Float32[comp.minor / (2 * sqrt(2 * log(2))) for comp in comps]
    # PA is measured from North
    pas = Float32[comp.pa + π/2 for comp in comps]

    if gpu
        @info "Sending model to GPU..."
        elapsed = Base.@elapsed begin
            # Send data to GPU device
            model_d = CuArrays.fill(ComplexF32(0), 4, length(lambdas), length(times))
            uvws_d = CuArray(uvws)
            lambdas_d = CuArray(lambdas)
            timesteps_d = CuArray(timesteps)
            lmns_d = CuArray(lmns)
            fluxes_d = CuArray(fluxes)
            lmns_d = CuArray(lmns)
            fluxes_d = CuArray(fluxes)
            majors_d = CuArray(majors)
            minors_d = CuArray(minors)
            pas_d = CuArray(pas)

            @cuda blocks=256 threads=256 shmem=sizeof(Float32)*length(lmns) + sizeof(Float32)*length(lambdas) predictgaussian!(model_d, uvws_d, lambdas_d, timesteps_d, lmns_d, fluxes_d, majors_d, minors_d, pas_d)
        end
        @info "Sent model to GPU for prediction elapsed $elapsed"

        return model_d
    else
        model = zeros(ComplexF32, 4, length(lambdas), length(times))

        elapsed = Base.@elapsed let
            predictgaussian!(model, uvws, lambdas, timesteps, lmns, fluxes, majors, minors, pas)
        end
        @info "Predicted model (CPU) elapsed $elapsed"

        return model
    end
end
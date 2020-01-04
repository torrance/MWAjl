function predict(
        uvws::Array{T, 2},
        times::Array{T, 1},
        lambdas::Array{T, 1},
        comps::Array{Component, 1},
        beam::Union{Beam, Nothing},
        pos0::Position,
    ) where {T <: AbstractFloat}

    freqs = 299792458 ./ lambdas

    # Construct an time index into each row
    elapsed = @elapsed begin
        unique_times = unique(times)
        timeidxs = indexin(times, unique_times)
    end
    @debug "Calculated time index elapsed $elapsed"

    # Precalculate l, m, n coordinates
    # TODO: This only needs to be done once
    elapsed = @elapsed begin
        lmns = zeros(3, length(comps))
        for (compidx, comp) in enumerate(comps)
            lmns[:, compidx] .= lmn(comp, pos0)
        end
    end
    @debug "Calculated lmns elapsed $elapsed"

    # Precalculate flux matrices for each source, for each frequency
    elapsed = @elapsed begin
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
        elapsed = @elapsed for (timeidx, time) in enumerate(unique_times)
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
            elapsed = @elapsed begin
                jones = reshape(
                    beamjones(beam, freq, alts, azs),
                    4, length(comps), length(unique_times)
                )
            end
            @debug "Retrieved beam Jones values, elapsed $elapsed"
            # apparent = J A J^H
            elapsed = @elapsed @views for timeidx in axes(jones, 3), compidx in axes(jones, 2)
                AxBH!(tmp, fluxes[:, compidx, chan, timeidx], jones[:, compidx, timeidx])
                AxB!(fluxes[:, compidx, chan, timeidx], jones[:, compidx, timeidx], tmp)
            end
            @debug "Applied beam correction to fluxes, elapsed $elapsed"
        end
    end

    # Allocate model array
    model = zeros(ComplexF32, 4, length(lambdas), length(times))
    elapsed = @elapsed predictionloop!(model, uvws, lambdas, timeidxs, lmns, fluxes)
    @debug "Predicted model elapsed $elapsed"

    return model
end

@inbounds @views function predictionloop!(model, uvws, lambdas, timeidxs, lmns, fluxes)
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
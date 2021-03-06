function consumer(ch, jones, converged, args)
    try
        while true
            elapsed = Base.@elapsed begin
                (data, model, weights, ants1, ants2, chanblock, timeblock) = take!(ch)
            end
            @info "Reading data (timeblock $timeblock chanblock $chanblock) elapsed $elapsed"
            # Initialize with most recent converged solution
            for i in chanblock:-1:1
                if converged[chanblock, timeblock]
                    jones[:, :, chanblock, timeblock] .= jones[:, :, i, timeblock]
                    # Sanitize any NaN present and set to Identity
                    for antid in axes(jones, 2), pol in axes(jones, 1)
                        if !isfinite(jones[pol, antid, chanblock, timeblock])
                            jones[pol, antid, chanblock, timeblock] = Int(pol in [1, 4])
                        end
                    end
                    break
                end
            end
            elapsed = Base.@elapsed converged[chanblock, timeblock], iterations = calibrate!(view(jones, :, :,  chanblock, timeblock), data, model, weights, ants1, ants2, args["max-iterations"], args["tolerance"]...)
            @info "Calibration (timeblock $timeblock chanblock $chanblock) complete, $iterations iterations, elapsed $elapsed"
        end
    catch e
        if isa(e, InvalidStateException)
            # Channel is closed, time to return our results
        else
            @error e
            rethrow(e)
        end
    end
    return jones
end

function producer(ch, mset, comps, beam, args)
    timeblocks = cld(mset.ntimesteps, args["timewidth"])
    chanblocks = cld(mset.nchans, args["chanwidth"])

    for timeblock in 1:timeblocks
        t1 = (timeblock - 1) * args["timewidth"] + 1
        t2 = min(mset.ntimesteps, timeblock * args["timewidth"])
        submset = taql("
            select * from \$1 where
            ANTENNA1 <> ANTENNA2
            and not FLAG_ROW
            and SUMSQR(UVW) > $(args["minuv"]^2)
            and SUMSQR(UVW) < $(args["maxuv"]^2)
            and TIME >= $(mset.unique_timesteps[t1])
            and TIME <= $(mset.unique_timesteps[t2])
        ", mset.table)

        # Load antenna IDs
        # Add one to satisfy one-indexed convention
        ants1 = column(submset, "ANTENNA1") .+ 1
        ants2 = column(submset, "ANTENNA2") .+ 1

        @info "$(length(ants1)) rows"

        # Additional columns needed for prediction
        local uvws, times
        if args["model"] !== nothing
            uvws = convert(Array{Float32}, column(submset, "UVW"))
            times = convert(Array{Float32}, column(submset, "TIME"))
        end

        chanblock = 1
        batchsize = args["chanwidth"] * args["nbatch"]
        for batchstart in 1:batchsize:mset.nchans
            batchend = min(mset.nchans, batchstart + batchsize - 1)

            # Fetch calibration data
            @debug "Fetching new batch of data"
            elapsed = Base.@elapsed begin
                if args["model"] === nothing
                    model = column(submset, args["modelcolumn"], blc=[1, batchstart], trc=[4, batchend])
                else
                    model = Threads.@spawn predict(uvws, times, mset.freqs[batchstart:batchend], comps, beam, mset.phasecenter, gpu=args["gpu"])
                end
                data = column(submset, args["datacolumn"], blc=[1, batchstart], trc=[4, batchend])
                flag = column(submset, "FLAG", blc=[1, batchstart], trc=[4, batchend])
                weights = column(submset, "WEIGHT_SPECTRUM", blc=[1, batchstart], trc=[4, batchend])
            end
            @info "Finished fetching new batch of data, elapsed $elapsed"

            # Model is built asynchronously, either by CPU threads or GPU
            # Wait for it to finish.
            if args["model"] !== nothing
                elapsed = Base.@elapsed model = Array(fetch(model))
                @info "Additional time waiting for model, elapsed $elapsed"
            end

            # Flag and sanitize data (eg. set NaN or Inf to 0)
            elapsed = Base.@elapsed sanitize!(data, model, flag)
            flag = nothing  # Allow GC of flag
            @debug "New data flagged and sanitised, elapsed $elapsed"

            # Send data to consumers
            elapsed = Base.@elapsed begin
                for chstart in 1:args["chanwidth"]:(batchend - batchstart + 1)
                    chend = min(mset.nchans, chstart + args["chanwidth"] - 1)
                    put!(ch, (
                        data[:, chstart:chend, :],
                        model[:, chstart:chend, :],
                        weights[:, chstart:chend, :],
                        ants1, ants2, chanblock, timeblock
                    ))
                    chanblock += 1
                end
            end
            @debug "Data set to consumers, elapsed $elapsed"
        end
    end
end
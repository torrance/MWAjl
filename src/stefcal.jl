using Statistics: mean
using LinearAlgebra: SingularException, mul!

#=
Assumption is that all autocorrelations are removed from data
data = [2x2, channels, rows]
model = [2x2, channels, rows]
=#
function calibrate!(jones::AbstractArray{Complex{Float64}, 2},
                    data::Array{Complex{T}, 3},
                    model::Array{Complex{T}, 3},
                    weights::Array{T, 3},
                    ants1::Array{S, 1},
                    ants2::Array{S, 1},
                    imax::Int,
                    amin::Float64,
                    amax::Float64,
                   )::Tuple{Bool, Int} where {T <: AbstractFloat, S <: Integer}

    nants = size(jones, 2)
    newjones = similar(jones)
    top = similar(jones)
    bot = similar(jones)
    failed = zeros(Bool, nants)
    distances = similar(jones, Float64)
    z = zeros(ComplexF64, 4)

    # Weight data
    weights = minimum(weights, dims=1)  # Set constant weight per polarization group
    data .*= weights
    model .*= weights

    iteration = 0
    while iteration < imax
        iteration += 1

        fill!(top, 0)
        fill!(bot, 0)
        calibrationloop(data, model, jones, ants1, ants2, top, bot, z)

        for antid in 1:nants
            if failed[antid]
                continue
            end

            try
                @views AdivB!(newjones[:, antid], top[:, antid], bot[:, antid])
            catch e
                if isa(e, SingularException)
                    # We set failed antennas to 0 rather than NaN to avoid special checks
                    # for NaN during the inner loop (conditionals are more expensive than just
                    # doing the extra work), and to allow the use of fastmath.
                    failed[antid] = true
                    jones[:, antid] .= 0
                    newjones[:, antid] .= 0
                else
                    rethrow(e)
                end
            end
        end

        # More than 4 antenna need to be present to get a good solution
        if nants - sum(failed) <= 4
            break
        end

        # On every even iteration, we test for convergence
        # and also set the new gain solution as the average of the last two,
        # as per Stefcal. This speeds up convergence.
        if iseven(iteration)
            distances[:, .~failed] .= abs2.(newjones[:, .~failed] - jones[:, .~failed])
            @debug "Iteration $iteration, failed antennae $(sum(failed)), distances" mean(distances[:, .~failed], dims=2)
            jones .= 0.5 * (jones + newjones)

            # Exit early if we reach stopping threshold
            distance = maximum(mean(distances[:, .~failed], dims=2))
            if distance < amax
                break
            end
        else
            # On odd iterations, we simply update the Jones matrix with the new one
            jones .= newjones
        end
    end

    # Set failed antennas to NaN. Also set zero solutions as failed,
    # which might arise if data for an antenna was zero, for some reason.
    jones[:, failed] .= NaN
    jones[:, (sum(jones, dims=1) .== 0)[:]] .= NaN

    # Exit criterion, in order of precedence
    distance = maximum(mean(distances[:, .~failed], dims=2))
    # First, if only 4 or fewer antennae remain, mark the solution as trash
    if nants - sum(failed) <= 4
        @info "Too many antenna solutions failed ($(sum(failed))) after $iteration iterations, setting solution block as failed"
        jones[:, :] .= NaN
        return false, iteration
    # Second, if we never reached the minimum threshold level, mark the entire solution as failed
    elseif distance > amin
        @info "Solution block failed to converge after $iteration iterations, setting as failed for all antennas (distance = $distance)"
        jones[:, :] .= NaN
        return false, iteration
    # Third, we exceeded the minimum threshold level, but not the maximum (ie. we didn't break early)
    elseif distance > amax
        @info "Solution block converged but did not meet amax threshold after $iteration iterations (distance = $distance)"
        return true, iteration
    # Finally, we exceeded the maximum threshold level and broke the iterations early
    else
        @info "Solution block converged after $iteration iterations (distance = $distance)"
        return true, iteration
    end
end

@inline @inbounds @views @fastmath function calibrationloop(data, model, jones, ants1, ants2, top, bot, z)
    for row in axes(data, 3)
        ant1 = ants1[row]
        ant2 = ants2[row]

       for chan in axes(data, 2)
            # Andre's calibrate: ( D J M^H ) / ( M J^H J M^H )
            AxBH!(z, jones[:, ant2], model[:, chan, row])  # J M^H
            plusAxB!(top[:, ant1], data[:, chan, row], z)  # D * z
            plusAHxB!(bot[:, ant1], z, z)

            AxB!(z, jones[:, ant1], model[:, chan, row])  # J (M^H)^H
            plusAHxB!(top[:, ant2], data[:, chan, row], z)  # D^H * z
            plusAHxB!(bot[:, ant2], z, z)

            # # Stefcal Paper: ( D^H J M ) / ( M^H J^H J M )
            # z = jones[:, :, ant2] * model[:, :, chan, row]'  # J M
            # top[:, :, ant1] += data[:, :, chan, row] * z  # D^H * z
            # bot[:, :, ant1] += z' * z  # z^H z

            # z = jones[:, :, ant1] * model[:, :, chan, row]  # J (M^H)
            # top[:, :, ant2] += data[:, :, chan, row]' * z  # (D^H)^H * z
            # bot[:, :, ant2] += z' * z  # z^H z
        end
    end
end

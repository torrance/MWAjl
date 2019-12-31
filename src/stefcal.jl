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
                   ) where {T <: AbstractFloat, S <: Integer}

    newjones = similar(jones)
    top = similar(jones)
    bot = similar(jones)
    failed = zeros(Bool, size(jones, 2))
    distances2 = similar(jones, Float64)
    z = zeros(ComplexF64, 4)
    amin2 = amin^2
    amax2 = amax^2
    
    for iteration in 1:imax
        fill!(top, 0)
        fill!(bot, 0)

        innerloop(data, model, jones, ants1, ants2, top, bot, z)

        for antid in axes(jones, 2)
            if failed[antid]
                continue
            end

            try
                @views AdivB!(newjones[:, antid], top[:, antid], bot[:, antid])
            catch e
                if isa(e, SingularException)
                    failed[antid] = true
                    jones[:, antid] .= 0
                    newjones[:, antid] .= 0
                else
                    rethrow(e)
                end
            end
        end

        # More than 4 antenna need to be present to get a good solution
        if sum(failed) + 4 >= size(jones, 2)
            @info "Solution block has too many failed antenna ($(sum(failed))) to continue, marking as failed after $iteration iterations"
            jones[:, :] .= NaN
            return
        end
    
        # On every even iteration, we test for convergence
        # and also set the new gain solution as the average of the last two,
        # as per Stefcal. This speeds up convergence.
        if iseven(iteration)
            distances2[:, .~failed] .= abs2.(newjones[:, .~failed] - jones[:, .~failed])
            jones .= 0.5 * (jones + newjones)

            # Exit early if we reach stopping threshold
            distance2 = maximum(mean(distances2[:, .~failed], dims=2))
            if distance2 < amax2
                @info "Solution block converged to amax threshold after $iteration iterations (distance = $(sqrt(distance2)))"
                break
            end
        else
            jones .= newjones
        end
    end

    # If the mean of the whole calibration block fails to meet amin, set it as failed
    distance2 = maximum(mean(distances2[:, .~failed], dims=2))
    if distance2 > amin2
        @info "Solution block failed to converge after $imax iterations, setting as failed for all antennas (distance = $(sqrt(distance2)))"
        jones[:, :] .= NaN
    elseif distance2 > amax2
        @info "Solution block converged but did not meet amax threshold after $imax iterations (distance = $(sqrt(distance2)))"
    end
    # Set any individual antennas that failed to meet the minimum threshold to NaN
    for (distance2, idx) in zip(findmax(distances2, dims=1)...)
        if distance2 > amin2
            antid = idx[2]
            @debug "Antenna $antid failed to converge" distance=sqrt(distance2)
            jones[:, antid] .= NaN
        end
    end
    # And finally set any failed antennas to NaN
    jones[:, failed] .= NaN
end

@inline @inbounds @views function innerloop(data, model, jones, ants1, ants2, top, bot, z)
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
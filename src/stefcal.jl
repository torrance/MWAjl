using Statistics: mean
using LinearAlgebra: SingularException

#=
Assumption is that all autocorrelations are removed from data
data = [2x2, channels, rows]
model = [2x2, channels, rows]
=#
function calibrate!(jones::AbstractArray{Complex{Float64}, 3},
                    data::Array{Complex{T}, 4},
                    model::Array{Complex{T}, 4},
                    weights::Array{T, 4},
                    ants1::Array{S, 1},
                    ants2::Array{S, 1},
                    amin::Float64,
                    amax::Float64,
                   ) where {T <: AbstractFloat, S <: Integer}

    newjones = similar(jones)
    top = similar(jones)
    bot = similar(jones)
    failed = zeros(Bool, size(jones, 3))
    distances = similar(jones, Float64)
    
    for iteration in 1:100
        fill!(top, 0)
        fill!(bot, 0)

        @views for row in axes(data, 4); 
            ant1 = ants1[row]
            ant2 = ants2[row]

            # TODO: Remove these conditionals
            # and set Jones = 0 for failed antennas?
            if ant1 == ant2 || failed[ant1] || failed[ant2]
                continue
            end

           for chan in axes(data, 3)
                # Andre's calibrate: ( D J M^H ) / ( M J^H J M^H )
                z = jones[:, :, ant2] * model[:, :, chan, row]'  # J M^H
                top[:, :, ant1] += data[:, :, chan, row] * z  # D * z
                bot[:, :, ant1] += z' * z

                z = jones[:, :, ant1] * model[:, :, chan, row]  # J (M^H)^H
                top[:, :, ant2] += data[:, :, chan, row]' * z  # D^H * z
                bot[:, :, ant2] += z' * z

                # # Stefcal Paper: ( D^H J M ) / ( M^H J^H J M )
                # z = jones[:, :, ant2] * model[:, :, chan, row]'  # J M
                # top[:, :, ant1] += data[:, :, chan, row] * z  # D^H * z
                # bot[:, :, ant1] += z' * z  # z^H z

                # z = jones[:, :, ant1] * model[:, :, chan, row]  # J (M^H)
                # top[:, :, ant2] += data[:, :, chan, row]' * z  # (D^H)^H * z
                # bot[:, :, ant2] += z' * z  # z^H z
            end
        end

        for antid in axes(jones, 3)
            if failed[antid]
                continue
            end

            try
                newjones[:, :, antid] = top[:, :, antid] / bot[:, :, antid]
            catch e
                if isa(e, SingularException)
                    failed[antid] = true
                else
                    rethrow(e)
                end
            end
        end

        # distances = abs.(newjones[:, :, .~failed] - jones[:, :, .~failed])
        # println("Max dist: ", maximum(mean(distances, dims=3)))
        # jones .= 0.25 * newjones + 0.75 * jones

        # if maximum(mean(distances, dims=3)) < 1E-7
        #     converged += 1
        #     if converged == 3
        #         break
        #     end
        # else
        #     convered = 0
        # end
    
        # On every even iteration, we test for convergence
        # and also set the new gain solution as the average of the last two,
        # as per Stefcal. This speeds up convergence.
        if iseven(iteration)
            distances[:, :, .~failed] .= abs.(newjones[:, :, .~failed] - jones[:, :, .~failed])
            jones .= 0.5 * (jones + newjones)

            # Exit early if we reach stopping threshold
            if maximum(mean(distances[:, :, .~failed], dims=3)) < amax
                @info "Solution converged to amax after $iteration iterations"
                break
            end
        else
            jones .= newjones
        end
    end

    # If we didn't meet the minimum tolerance level, set the whole solution to NaN
    # Otherwise, just set those specific entries that failed to NaN
    if maximum(mean(distances[:, :, .~failed], dims=3)) > amin
        @info "Solution did not converge to amin" pdistances=mean(distances[:, :, .~failed], dims=3)
        jones .= NaN
    else
        jones[:, :, failed] .= NaN
    end

    return jones
end
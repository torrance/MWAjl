using Statistics: mean
using LinearAlgebra: SingularException, mul!

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
                    imax::Int,
                    amin::Float64,
                    amax::Float64,
                   ) where {T <: AbstractFloat, S <: Integer}

    newjones = similar(jones)
    top = similar(jones)
    bot = similar(jones)
    failed = zeros(Bool, size(jones, 3))
    distances2 = similar(jones, Float64)
    z = zeros(ComplexF64, 2, 2)
    tmp = zeros(ComplexF64, 2, 2)
    amin2 = amin^2
    amax2 = amax^2
    
    for iteration in 1:imax
        fill!(top, 0)
        fill!(bot, 0)

        innerloop(data, model, jones, ants1, ants2, top, bot, z, tmp)

        for antid in axes(jones, 3)
            if failed[antid]
                continue
            end

            try
                newjones[:, :, antid] .= top[:, :, antid] / bot[:, :, antid]
            catch e
                if isa(e, SingularException)
                    failed[antid] = true
                    jones[:, :, antid] .= 0
                else
                    rethrow(e)
                end
            end
        end
    
        # On every even iteration, we test for convergence
        # and also set the new gain solution as the average of the last two,
        # as per Stefcal. This speeds up convergence.
        if iseven(iteration)
            distances2[:, :, .~failed] .= abs2.(newjones[:, :, .~failed] - jones[:, :, .~failed])
            jones .= 0.5 * (jones + newjones)

            # Exit early if we reach stopping threshold
            distance2 = maximum(mean(distances2[:, :, .~failed], dims=3))
            if distance2 < amax2
                @info "Solution block converged to amax threshold after $iteration iterations (distance = $(sqrt(distance2)))"
                break
            end
        else
            jones .= newjones
        end
    end

    # If the mean of the whole calibration block fails to meet amin, set it as failed
    distance2 = maximum(mean(distances2[:, :, .~failed], dims=3))
    if distance2 > amin2
        @info "Solution block failed to converge after $imax iterations, setting as failed for all antennas (distance = $(sqrt(distance2)))"
        jones[:, :, :] .= NaN
    elseif distance2 > amax2
        @info "Solution block converged but did not meet amax threshold after $imax iterations (distance = $(sqrt(distance2)))"
    end
    # Set any individual antennas that failed to meet the minimum threshold to NaN
    for (distance2, idx) in zip(findmax(distances2, dims=[1, 2])...)
        if distance2 > amin2
            antid = idx[3]
            @debug "Antenna $antid failed to converge" distance=sqrt(distance2)
            jones[:, :, antid] .= NaN
        end
    end
    # And finally set any failed antennas to NaN
    jones[:, :, failed] .= NaN
end

@inline @inbounds @views function innerloop(data, model, jones, ants1, ants2, top, bot, z, tmp)
    for row in axes(data, 4)
        ant1 = ants1[row]
        ant2 = ants2[row]

       for chan in axes(data, 3)
            # Andre's calibrate: ( D J M^H ) / ( M J^H J M^H )
            mul!(z, jones[:, :, ant2], model[:, :, chan, row]')  # J M^H
            mul!(tmp, data[:, :, chan, row], z)
            top[:, :, ant1] .+=  tmp # D * z
            mul!(tmp, z', z)
            bot[:, :, ant1] .+= tmp

            mul!(z, jones[:, :, ant1], model[:, :, chan, row])  # J (M^H)^H
            mul!(tmp, data[:, :, chan, row]', z)
            top[:, :, ant2] .+= tmp  # D^H * z
            mul!(tmp, z', z)
            bot[:, :, ant2] .+= tmp

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
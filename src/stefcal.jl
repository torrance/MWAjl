using Statistics: mean
using LinearAlgebra: SingularException

#=
data = [2x2, channels, rows]
model = [2x2, channels, rows]
=#
function calibrate(data::Array{Complex{T}, 4}, model::Array{Complex{T}, 4}, ants1::Array{S, 1}, ants2::Array{S, 1}) where {T <: AbstractFloat, S <: Integer}
    nants =  max(maximum(ants1), maximum(ants2))
    jones = zeros(Complex{T}, 2, 2, nants)
    newjones = similar(jones)
    top = similar(jones)
    bot = similar(jones)
    failed = zeros(Bool, nants)

    # Initialize Jones matrices to Identity
    jones[1, 1, :] .= 1
    jones[2, 2, :] .= 1

    print("Solving...")
    for iteration in 1:100
        print("*")
        fill!(top, 0)
        fill!(bot, 0)

        for row in axes(data, 4), chan in axes(data, 3)
            ant1 = ants1[row]
            ant2 = ants2[row]

            if ant1 == ant2 || failed[ant1] || failed[ant2]
                continue
            end

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
        print("\b.")

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
            distances = abs.(newjones[:, :, .~failed] - jones[:, :, .~failed])
            jones .= 0.5 * (jones + newjones)

            # Exit early if we reach stopping threshold
            if maximum(mean(distances, dims=3)) < 1E-7
                break
            end
        else
            jones .= newjones
        end
    end
    print("Done\n")

    jones[:, :, failed] .= NaN
    return jones
end
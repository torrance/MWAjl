
function sanitize!(data::Array{Complex{T}, 3}, model::Array{Complex{T}, 3}, flag::Array{Bool, 3}) where T <: AbstractFloat
    for row in axes(data, 3), chan in axes(data, 2)
        for pol in axes(data, 1)
            if flag[pol, chan, row] || !isfinite(data[pol, chan, row]) || !isfinite(model[pol, chan, row])
                data[:, chan, row] .= 0
                model[:, chan, row] .= 0
                break  # We've flagged all polarizations, no need to check the others
            end
        end
    end
end
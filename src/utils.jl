
function sanitize!(data::Array{Complex{T}, 4}, model::Array{Complex{T}, 4}, flag::Array{Bool, 4}) where T <: AbstractFloat
    for row in axes(data, 4), chan in axes(data, 3)
        for polx in axes(data, 2), poly in axes(data, 1)
            if flag[polx, poly, chan, row] || !isfinite(data[polx, poly, chan, row]) || !isfinite(model[polx, poly, chan, row])
                data[:, :, chan, row] .= 0
                model[:, :, chan, row] .= 0
                break  # We've flagged all polarizations, no need to check the others
            end
        end
    end
end
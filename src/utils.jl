
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

function polyfit(x::Array{Float64}, y::Array{Float64}, n::Int)
    length(x) == length(y) || throw(DomainError)
    1 <= n <= length(x) - 1 || throw(DomainError)

    A = Array{Float64}(undef, length(x), n+1)
    A[:, 1] .= 1
    for i in 1:n
        A[:, i+1] .= A[:, i] .* x
    end
    A \ y
end
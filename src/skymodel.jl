struct Position
    ra::Float64
    dec::Float64
end

struct Measurement
    ν::Float64
    stokes::Array{Float64}  # [I Q U V]
end

abstract type Spectrum end

struct Component
    position::Position
    type::String
    spectrum::Spectrum
end

# TODO: Ideally, want a type alias for Array{Measurements}
struct Measurements <: Spectrum
    ms::Array{Measurement}
end

struct SED <: Spectrum
    ν::Float64
    stokes::Array{Float64}  # I Q U V
    coeffs::Array{Float64}
end

struct Source
    name::String
    components::Array{Component}
end


function instrumental(comp::Component)::Array{Float64}
    I, Q, U, V = stokes(comp.spectrum, ν)

    # Calculate instrumental polarization
    XX = I + Q
    XY = U + 1im * V
    YX = U - 1im * V
    YY = I - Q

    # TODO: Add beam
    return [XX XY; YX YY]
end


function stokes(sed::SED, ν::Float64)::Array{Float64}
    # Calculate the frequency adjust Stokes parameters
    stokes = zeros(4)  # [ I Q U V ]
    for (i, stoke) in enumerate(sed.stokes)
        if stoke <= 0
            stokes[i] = 0
            continue
        end

        # Flux at new frequency is given by:
        # logS = logS_0 + alpha * ( log nu - log nu_0) + beta * (log^2 nu - log^2 nu_0) ..
        logS = log(stoke)
        for (power, coeff) in enumerate(sed.coeffs)
            logS += coeff * (log(ν)^power - log(sed.ν)^power)
        end
        stokes[i] =  ℯ^logS
    end

    stokes
end


function stokes(measurements::Measurements, ν::Float64)::Array{Float64}
    # We interpolate in two ways
    # 1. If freq lies between two measured values, we linearly interpolate (in log spaace)
    #    between these values.
    # 2. Otherwise, we calculate spectral index using all measured values.

    stokes = zeros(4)  # [I Q U V]
    ms = measurements.ms

    # Attempt Case 1
    for i in 1:length(ms)-1
        if ms[i].ν <= ν < ms[i + 1].ν
            # We've found a pair of measurements adjacent to our requested frequency.
            # Loop over stokes parameters
            # TODO: Handle case where one of the measurements is <= 0
            for j in 1:4
                println(ms[i].stokes[j], " ", ms[i + 1].stokes[j])
                if ms[i].stokes[j] ≤ 0 && ms[i + 1].stokes[j] ≤ 0
                    stokes[j] = 0
                else
                    c, α = polyfit(log.([ms[i].ν, ms[i + 1].ν]), log.([ms[i].stokes[j], ms[i + 1].stokes[j]]), 1)
                    stokes[j] = ℯ^(c + α * log(ν))
                end
            end
            return stokes
        end
    end

    # Fallback: Case 2
    # Loop over stokes paramters
    for i in 1:4
        logνs = Float64[]
        logSs = Float64[]
        for m in ms
            if m.stokes[i] > 0
                push!(logνs, log(m.ν))
                push!(logSs, log(m.stokes[i]))
            end
        end

        # Handle case if all measurements are zero
        if length(logSs) == 0
            stokes[i] = 0
        else
            # TODO: handle case when we only have 1 measurement
            c, α = polyfit(logνs, logSs, 1)
            stokes[i] = ℯ^(c + α * log(ν))
        end
    end

    stokes
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
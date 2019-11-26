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

struct Measurements <: Spectrum
    ms::Array{Measurement}
    function Measurements(ms::Array{Measurement})
        if length(ms) == 0
            throw(ArgumentError("Require at least one measurement"))
        end
        ms = sort(ms, by=m -> m.ν)
        new(ms)
    end
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


function Base.getindex(ms::Measurements, i)
    ms.ms[i]
end


function Base.length(ms::Measurements)
    length(ms.ms)
end


function Base.iterate(ms::Measurements)
    iterate(ms.ms)
end


function Base.iterate(ms::Measurements, state)
    iterate(ms.ms, state)
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


function stokes(ms::Measurements, ν::Float64)::Array{Float64}
    # We interpolate in three ways
    # 1. If just one value is provide, we assume 0 spectral index (ie. constant)
    # 2. If freq lies between two measured values, we linearly interpolate (in log spaace)
    #    between these values.
    # 3. Otherwise, we calculate spectral index using all measured values.

    stokes = zeros(4)  # [I Q U V]

    # Case 1
    if length(ms) == 1
        return ms[1].stokes[:]
    end

    # Attempt case 2
    for i in 1:length(ms)-1
        if ms[i].ν <= ν <= ms[i + 1].ν
            # We've found a pair of measurements adjacent to our requested frequency.
            # Loop over stokes parameters
            # TODO: Handle case where one of the measurements is <= 0
            for j in 1:4
                # If both measurements are 0, we set to 0
                if ms[i].stokes[j] ≤ 0 && ms[i + 1].stokes[j] ≤ 0
                    stokes[j] = 0
                # If one of the measurements is 0, then log space interpolation will fail
                # so we fallback to linear space interpolation
                elseif ms[i].stokes[j] ≤ 0 || ms[i + 1].stokes[j] ≤ 0
                    c, m = polyfit([ms[i].ν, ms[i + 1].ν], [ms[i].stokes[j], ms[i + 1].stokes[j]], 1)
                    stokes[j] = c + m * ν
                # Interpolate in log space
                else
                    c, α = polyfit(log.([ms[i].ν, ms[i + 1].ν]), log.([ms[i].stokes[j], ms[i + 1].stokes[j]]), 1)
                    stokes[j] = ℯ^(c + α * log(ν))
                end
            end
            return stokes
        end
    end

    # Fallback: Case 3
    # Loop over stokes parameters
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
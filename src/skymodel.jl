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


function instrumental(comp::Component, ν::T)::Array{T} where T <: AbstractFloat
    I, Q, U, V = stokes(comp.spectrum, ν)

    # Calculate instrumental polarization
    XX = I + Q
    XY = U + 1im * V
    YX = U - 1im * V
    YY = I - Q

    # TODO: Add beam
    return [XX, XY, YX, YY]
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


"""
    stokes(ms::Measurements, ν::Float64)

Given a Measurements object (containing one more more Measurement objects),
interpolate each of Stokes I, Q, U, V to the frequency ν.

This function has been implemented to identically match the behaviour of the original `calibrate`.

We interpolate in three ways:

1. If just one measurement is provide, we assume a flat spectral index (i.e. constant).
2. If the frequency ν lies between two measured values, we linearly interpolate (in log/log space)
   between these values. EXCEPT if either of the two neighbour values is <= 0, then we
   linearly interpolate in linear space between the two neighbouring measurements.
3. Otherwise, we linearly interpolate in log/log space based on a least squares fit to _all_
    strictly positive measurements. If none are strictly positive, we return 0. If just one is strictly positive, the spectrum is assumed flat.

!!! warning
    There are plenty of pathological combinations of measurement values that will give unexpected results. For example, given measurements 100 MHz => 1 Jy and 200 MHz => 0 Jy, then requesting the flux density at 300 MHz will give 1 Jy.

    To avoid such unexpected results, ensure all measurements are strictly positive (i.e. > 0).

    Or for complete control over interpolation, it is recommended to instead provide an SED.
"""
function stokes(ms::Measurements, ν::Float64)::Array{Float64}
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
            for j in 1:4
                # If one of the measurements is <= 0, then log space interpolation will fail
                # so we fallback to linear space interpolation
                if ms[i].stokes[j] ≤ 0 || ms[i + 1].stokes[j] ≤ 0
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
        elseif length(logSs) == 1
            stokes[i] = exp(logSs[1])
        else
            c, α = polyfit(logνs, logSs, 1)
            stokes[i] = ℯ^(c + α * log(ν))
        end
    end

    stokes
end


function lmn(comp::Component, pos0)
    ra, dec = comp.position.ra, comp.position.dec
    ra0, dec0 = pos0.ra, pos0.dec

    l = cos(dec) * sin(ra - ra0)
    m = sin(dec) * cos(dec0) - cos(dec) * sin(dec0) * cos(ra - ra0)
    n = sqrt(1 - l^2 - m^2)
    return l, m, n
end
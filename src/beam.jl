using Pkg.Artifacts

const mwabeam = joinpath(artifact"mwabeam", "mwa_full_embedded_element_pattern.h5")
const hyperbeam = "libmwa_hyperbeam"

mutable struct Beam
    ptr::Ptr
    delays::Array{Int32, 1}
    amps::Array{Float64, 1}
end

function Beam(delays::Array{Int32, 1}, amps=nothing::Union{Nothing, Array{Float64, 1}})
    if amps === nothing
        amps = ones(size(delays))
    end
    beam_ptr = @ccall hyperbeam.new_fee_beam(mwabeam::Cstring)::Ptr
    beam = Beam(beam_ptr, delays, amps)
    finalizer(beam) do f
        @ccall hyperbeam.free_fee_beam(beam.ptr::Ptr)::Cvoid
        return beam
    end
    return beam
end

"""
    calc_jones()

    az: azimuth, measured North through East [radian]
    za: zenith angle [radian]
    freq: frequency [Hz]
"""
function calc_jones(beam::Beam, az::Float64, za::Float64, freq::Int, zenithnorm=true)
    jones_ptr = @ccall hyperbeam.calc_jones(beam.ptr::Ptr, az::Cdouble, za::Cdouble, freq::Cuint, beam.delays::Ptr{Cuint}, beam.amps::Ptr{Cdouble}, zenithnorm::Cuchar)::Ptr{Cdouble}
    jones = unsafe_wrap(Array, jones_ptr, 8, own=true)

    # Convert pairs of doubles to complex
    return jones[1:2:end] .+ 1im * jones[2:2:end]
end

"""
    calc_jones()

    azs: azimuth, measured North through East [radian] Array[N]
    zas: zenith angle [radian] Array[N]
    freq: frequency [Hz]

    Returns: Array[4, N]
"""
function calc_jones(beam::Beam, azs::Array{Float64}, zas::Array{Float64}, freq::UInt32, zenithnorm=true)
    @assert length(azs) == length(zas)
    N = length(azs)
    jones_ptr = @ccall hyperbeam.calc_jones_array(beam.ptr::Ptr, N::Cuint, azs::Ptr{Cdouble}, zas::Ptr{Cdouble}, freq::Cuint, beam.delays::Ptr{Cuint}, beam.amps::Ptr{Cdouble}, zenithnorm::Cuchar)::Ptr{Cdouble}
    jones = unsafe_wrap(Array, jones_ptr, N * 8, own=true)
    jones = reshape(jones, 8, N)

    # Convert pairs of doubles to complex
    return jones[1:2:end, :] .+ 1im * jones[2:2:end, :]
end

function closest_freq(beam::Beam, freq::UInt32)
    return @ccall hyperbeam.closest_freq(beam.ptr::Ptr, freq::Cuint)::Cuint
end
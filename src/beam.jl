using Pkg.Artifacts

mutable struct Beam
    delays::Array{Int32, 1}
    path::String
    ptr::Ptr{Cvoid}
end

const libbeam = joinpath(artifact"libbeam.so", "libbeam.so")

function Beam(delays::Array{Int32, 1}, path::String)
    amps = ones(16)
    ptr = ccall((:beam_new, libbeam), Ptr{Cvoid}, (Ptr{Cdouble}, Ptr{Cdouble}, Cstring), convert(Array{Cdouble}, delays), amps, path)
    beam = Beam(delays, path, ptr)
    finalizer(del, beam)
    return beam
end

function del(beam::Beam)
    ccall((:beam_del, libbeam), Cvoid, (Ptr{Cvoid},), beam.ptr)
end

function beamjones(beam::Beam, freq::Float64, alts::Array{Float64, 1}, azs::Array{Float64, 1})
    if length(alts) != length(azs)
        throw(ArgumentError("alts and azs must be the same length"))
    end
    N = length(alts)

    jones = Array{ComplexF64, 2}(undef, 4, N)
    ccall((:beamjones, libbeam), Cvoid, (Ptr{Cvoid}, Cint, Csize_t, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{ComplexF64}), beam.ptr, round(Cint, freq), N, alts, azs, jones)
    return jones
end

function closest_freq(beam::Beam, freq::Float64)
    closest = ccall((:find_closest_freq, libbeam), Cint, (Ptr{Cvoid}, Cint), beam.ptr, round(Cint, freq))
    return convert(Float64, closest)
end
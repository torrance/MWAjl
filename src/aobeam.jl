mutable struct AOBeam
    delays::Array{Int32, 1}
    path::String
    ptr::Ptr{Cvoid}
end

const libaobeam = string(@__DIR__, "/libaobeam.so")

function AOBeam(delays::Array{Int32, 1}, path::String)
    amps = ones(16)
    ptr = @threadcall((:beam_new, libaobeam), Ptr{Cvoid}, (Ptr{Cdouble}, Ptr{Cdouble}, Cstring), convert(Array{Cdouble}, delays), amps, path)
    beam = AOBeam(delays, path, ptr)
    finalizer(del, beam)
    return beam
end

function del(beam::AOBeam)
    ccall((:beam_del, libaobeam), Cvoid, (Ptr{Cvoid},), beam.ptr)
end

function beamjones(beam::AOBeam, freq::Float64, alts::Array{Float64, 1}, azs::Array{Float64, 1})
    if length(alts) != length(azs)
        throw(ArgumentError("alts and azs must be the same length"))
    end
    N = length(alts)

    jones = Array{ComplexF64, 2}(undef, 4, N)
    @threadcall((:beamjones, libaobeam), Cvoid, (Ptr{Cvoid}, Cint, Csize_t, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{ComplexF64}), beam.ptr, round(Cint, freq), N, alts, azs, jones)
    return jones
end

function closest_freq(beam::AOBeam, freq::Float64)
    closest = @threadcall((:find_closest_freq, libaobeam), Cint, (Ptr{Cvoid}, Cint), beam.ptr, round(Cint, freq))
    return convert(Float64, closest)
end
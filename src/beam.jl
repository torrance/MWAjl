# For the meantime, this is just a wrapper of Python mwabeam
using PyCall

struct Beam
    delays::Array{Int32, 1}
end

function beamjones(beam::Beam, freq, alt, az)
    jones = mwapb.MWA_Tile_full_EE(π/2 .- alt, az, freq, delays=beam.delays, jones=true, interp=false)
    # The output of mwa_pb is [n, 2, 2], where n is length(alt) = length(az)
    # However, we want this in [4, n] format. Note also we swap the indices of the
    # 2x2 matrices to account for the fact that Python is row major and Julia is column major.
    jones = permutedims(jones, (3, 2, 1))
    jones = reshape(jones, 4, size(jones, 3))
    return jones
end
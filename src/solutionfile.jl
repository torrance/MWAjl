# The solution file is used for storing calibration Jones matrices. The format is as follows:
#  Bytes |  Description
# -------+---------------
#  0- 7  |  string intro ; 8-byte null terminated string "MWAOCAL"
#  8-11  |  int fileType ; always 0, reserved for indicating something other than complex Jones solutions
# 12-15  |  int structureType ; always 0, reserved for indicating different ordering
# 16-19  |  int intervalCount ; Number of solution intervals in file
# 20-23  |  int antennaCount ; Number of antennas that were in the measurement set (but were not necessary all solved for)
# 24-27  |  int channelCount ; Number of channels in the measurement set
# 28-31  |  int polarizationCount ; Number of polarizations solved for -- always four.
# 32-39  |  double startTime ; Start time of solutions (AIPS time)
# 40-47  |  double endTime ; End time of solutions (AIPS time)
# -------+-------------------
# After the header follow 2 x nSolution doubles, with

# nSolutions = nIntervals * nAntennas * nChannels * nPols

# Ordered in the way as given, so:
# double 0 : real of interval 0, antenna 0, channel 0, pol 0
# double 1 : imaginary of interval 0, antenna 0, channel 0, pol 0
# double 2 : real of interval 0, antenna 0, channel 0, pol 1
# ...
# double 8 : real of interval 0, antenna 0, channel 1, pol 0
# double nChannel x 8 : real of interval 0, antenna 1, channel 0, pol 0
# etc.

# here, ints are always 32 bits unsigned integers, doubles are IEEE double precision 64 bit floating points.
# If a solution is not available, either because no data were selected during calibration for this interval
# or because the calibration diverged, a "NaN" will be stored in the doubles belonging to that solution.
function writesolution(io::IO, solution::Array{ComplexF64, 4}, start::Float64, finish::Float64)
    pol, nants, nchans, intervals = size(solution)
    solution = permutedims(solution, [1, 3, 2, 4])

    write(io, UInt8['M', 'W', 'A', 'O', 'C', 'A', 'L', 0])  # intro
    write(io, reinterpret(UInt8, Int32[0]))                 # fileType
    write(io, reinterpret(UInt8, Int32[0]))                 # structureType
    write(io, reinterpret(UInt8, Int32[intervals]))         # intervalCount
    write(io, reinterpret(UInt8, Int32[nants]))             # antennaCount
    write(io, reinterpret(UInt8, Int32[nchans]))            # channelCount
    write(io, reinterpret(UInt8, Int32[pol]))               # polarizationCount
    write(io, reinterpret(UInt8, Float64[start]))           # startTime
    write(io, reinterpret(UInt8, Float64[finish]))          # endTime
    write(io, reinterpret(UInt8, solution))
end
using ArgParse
using Dates
using Distributed
using FileIO
using JLD
using Logging

s = ArgParseSettings()
@add_arg_table s begin
    "--model", "-m"
        help="The path to the model file (aoskymodel 1.2 format). If absent, will use MODEL_DATA column (ie. self-calibration)."
    "--chanwidth", "-c"
        help="Find calibration solutions for blocks of channels of `chanwidth`. Set to 0 to find a single solution for all channels."
        arg_type=Int
        default=1
    "--timewidth", "-t"
        help="Find calibrations solutions for blocks of time steps of `timewidth`. Defaults to 0, which implies an infinite time width."
        arg_type=Int
        default=0
    "--minuv"
        help="Exclude baselines shorter than this length from affecting calibration solution (metres)."
        arg_type=Float64
        default=0.0
    "--maxuv"
        help="Exclude baselines longer than this length from affecting calibration solutions (metres)."
        arg_type=Float64
        default=9E99
    "--tolerance", "-a"
        help="These two values determine whether 1) a solution has sufficiently converged that we can accept its answer and 2) whether it has converged enough that we may stop prior to reaching `max-iterations`. These values are tested after each iteration by comparing the magnitude difference between the current solution and the solution from the last iteration."
        arg_type=Float64
        nargs=2
        default=[1E-5, 1E-8]
    "--max-iterations", "-i"
        help="The maximum number of iterations allowed when solving for a single solution unit (ie. for a given channel and time block)."
        arg_type=Int
        default=50
    "--apply-beam"
        help="The path to the MWA HD5 beam file."
        arg_type=String
    "--datacolumn"
        arg_type=String
        default="DATA"
    "--modelcolumn"
        arg_type=String
        default="MODEL_DATA"
    "--nbatch"
        help="When reading from `datacolumn` and `modelcolumn`, this parameter specifies how many `chanwidth` arrays to read."
        arg_type=Int
        default=1
    "--processes", "-j"
        help="Maximum simulaneous CPU processes to use. Default value (-1) will set to number of logical CPU cores. 0 will run synchronously."
        arg_type=Int
        default=-1
    "--verbose"
        action=:store_true
    "--debug"
        action=:store_true
    "mset"
        help="The path to the measurement set to be calibrated."
        required=true
    "solution"
        help="The output solution file name and path. eg. folder/to/go/solution.bin"
        required=true
end
args = parse_args(s)

# Set logging levels
if args["debug"]
    global_logger(Logging.ConsoleLogger(stdout, Logging.Debug))
elseif args["verbose"]
    global_logger(Logging.ConsoleLogger(stdout, Logging.Info))
else
    global_logger(Logging.ConsoleLogger(stdout, Logging.Warn))
end

# Create workers and load MWAjl
if args["processes"] < 0
    @info "Defaulting to use $(Sys.CPU_THREADS) processes"
    args["processes"] = Sys.CPU_THREADS
end
addprocs(args["processes"])
@everywhere using MWAjl

@info "Measurement set: $(args["mset"])"

# Lets get the metadata about this measurement set
# Frequency
spw = Table(args["mset"] * "/SPECTRAL_WINDOW")
freqs = column(spw, "CHAN_FREQ")
spw = nothing
if size(freqs, 2) > 1
    println("Fatal: There is more than one spectral window present in the measurement set. We only know how to handle one.")
    exit(1)
end
freqs = freqs[:, 1]
λs = 299792458 ./ freqs
nchans = length(freqs)
@info " Total channels: $nchans"

# Timesteps
mset = Table(args["mset"])
timesteps = sort(unique(column(mset, "TIME")))
ntimesteps = length(timesteps)
@info " Total timesteps: $ntimesteps"

# Antennas
antennas = Table(args["mset"] * "/ANTENNA")
nants = size(column(antennas, "POSITION"), 2)
antennas = nothing
@info " Number of antennas: $nants"

if args["chanwidth"] < 1
    args["chanwidth"] = length(freqs)
end
if length(freqs) % args["chanwidth"] != 0
    @warn "--chanwidth does not evenly divide the total number of channels"
end

if args["timewidth"] < 1
    args["timewidth"] = length(timesteps)
end
if length(timesteps) % args["timewidth"] != 0
    @warn "--timewidth does not evenly divide the total number of time steps"
end

# Initialize final Jones matrix to Identity
# [pol, pol, antid, channels, timesteps]
timeblocks = cld(ntimesteps, args["timewidth"])
chanblocks = cld(nchans, args["chanwidth"])
jones = zeros(Complex{Float64}, 2, 2, nants, chanblocks, timeblocks)
jones[1, 1, :, :, :] .= 1
jones[2, 2, :, :, :] .= 1
@debug "Created Jones array: ", size(jones)

# Initialise workers
ch = RemoteChannel(() -> Channel{
    Tuple{Array{ComplexF32, 4}, Array{ComplexF32, 4}, Vector{Int}, Vector{Int}, Int, Int}
}(args["nbatch"]))
futures = Future[]
for pid in workers()
    f = remotecall(function (ch, jones, args)
        try
            while true
                elapsed = @elapsed begin
                    (data, model, ants1, ants2, chanblock, timeblock) = take!(ch)
                end
                @info "Reading data (timeblock $timeblock chanblock $chanblock) elapsed $elapsed"
                elapsed = @elapsed calibrate!(view(jones, :, :, :,  chanblock, timeblock), data, model, similar(data, Float32), ants1, ants2, args["max-iterations"], args["tolerance"]...)
                @info "Calibration (timeblock $timeblock chanblock $chanblock) elapsed $elapsed"
            end
        catch e
            if isa(e, InvalidStateException)
                # Channel is closed, time to return our results
            else
                # TODO: propagate exception to main thread
            end
        end
        return jones
    end, pid, ch, jones, args)
    push!(futures, f)
end


for timeblock in 1:timeblocks
    t1 = (timeblock - 1) * args["timewidth"] + 1
    t2 = min(ntimesteps, timeblock * args["timewidth"])
    submset = taql("
        select * from \$1 where 
        ANTENNA1 <> ANTENNA2
        and not FLAG_ROW
        and SUM(SUMSQR(UVW)) > $(args["minuv"])
        and SUM(SUMSQR(UVW)) < $(args["maxuv"])
        and TIME >= $(timesteps[t1])
        and TIME <= $(timesteps[t2])
    ", mset)

    # Load antenna IDs
    # Add one to satisfy one-indexed convention
    ants1 = column(submset, "ANTENNA1") .+ 1
    ants2 = column(submset, "ANTENNA2") .+ 1

    chanblock = 1
    batchsize = args["chanwidth"] * args["nbatch"]
    for batchstart in 1:batchsize:nchans
        batchend = min(nchans, batchstart + batchsize - 1)

        # Read in data from measurement set
        @debug "Fetching new batch of data"
        elapsed = @elapsed begin
            data = column(submset, "DATA", blc=[1, batchstart], trc=[4, batchend])
            model = column(submset, "MODEL_DATA", blc=[1, batchstart], trc=[4, batchend])
            flag = column(submset, "FLAG", blc=[1, batchstart], trc=[4, batchend])
        end
        @debug "Finished fetching new batch of data, elapsed $elapsed"

        # Reshape so that Jones matrices from (4,) to (2, 2), as expected by calibrate
        data = reshape(data, 2, 2, size(data)[2:end]...)
        model = reshape(model, 2, 2, size(model)[2:end]...)
        flag = reshape(flag, 2, 2, size(flag)[2:end]...)

        # Flag and sanitize data (eg. set NaN or Inf to 0)
        elapsed = @elapsed sanitize!(data, model, flag)
        flag = nothing  # Allow GC of flag
        @debug "New data flagged and sanitised, elapsed $elapsed"

        # Send data to workers
        for chstart in 1:args["chanwidth"]:(batchend - batchstart + 1)
            chend = chstart + args["chanwidth"] - 1
            put!(ch, (data[:, :, chstart:chend, :], model[:, :, chstart:chend, :], ants1, ants2, chanblock, timeblock))
            chanblock += 1
        end
        data = nothing
        model = nothing
    end
end

# Wait on workers to complete and combine final jones matrix
close(ch)  # When the channel is drained, this signals to workers to return
for f in futures
    jones .+= fetch(f)
end

# Invert jones Matrices so that they can be applied as J D J^H
for timeblock in axes(jones[5])
    for chanblock in axes(jones[4])
        for antid in axes(jones[3])
            jones[:, :, antid, chanblock, timeblock] .= inv(jones[:, :, antid, chanblock, timeblock])
        end
    end
end

save(File(format"JLD", args["solution"]), "jones", jones)

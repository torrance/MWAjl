using ArgParse
using FileIO
using JLD
using Logging
using MWAjl

s = ArgParseSettings()
@add_arg_table s begin
    "--model", "-m"
        help="The path to the model file (aoskymodel 1.2 format). If absent, will use MODEL_DATA column (ie. self-calibration)."
        arg_type=String
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
    "--gpu"
        help="Use GPU acceleration for some calculations."
        action=:store_true
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
    global_logger(Logging.ConsoleLogger(stderr, Logging.Debug))
elseif args["verbose"]
    global_logger(Logging.ConsoleLogger(stderr, Logging.Info))
else
    global_logger(Logging.ConsoleLogger(stderr, Logging.Warn))
end

# Determine if we are running selfcalibration or predicting our own model
# If the latter, read in the model and construct a list of components
if args["model"] === nothing
    global comps = nothing
    @info "Self-calibration mode (using $(args["modelcolumn"]) column as model)"
else
    global comps = Component[]
    open(string(args["model"])) do f
        sources = parse_model(f)
        for source in sources, comp in source.components
            push!(comps, comp)
        end
    end
    @info "Sky model mode with $(length(comps)) components"
end

# Open the measurement set and load required metadata
mset = MeasurementSet(args["mset"])

# Load beam object if --apply-beam is set
if args["apply-beam"] !== nothing && mset.mwadelays !== nothing
    global beam = AOBeam(mset.mwadelays, args["apply-beam"])
elseif args["apply-beam"] !== nothing
    @error "Unable to load MWA tile delays from measurement set, which is required with --apply-beam"
    exit(1)
else
    global beam = nothing
end

# Set default channel width if --chanwidth==0
if args["chanwidth"] < 1 || args["chanwidth"] > mset.nchans
    args["chanwidth"] = mset.nchans
end
if mset.nchans % args["chanwidth"] != 0
    @warn "--chanwidth does not evenly divide the total number of channels"
end

# Set default time width if --timewidth==0
if args["timewidth"] < 1 || args["timewidth"] > mset.ntimesteps
    args["timewidth"] = mset.ntimesteps
end
if mset.ntimesteps % args["timewidth"] != 0
    @warn "--timewidth does not evenly divide the total number of time steps"
end

# Based on time width and chan width, calculate the number of timeblocks and chanblocks
timeblocks = cld(mset.ntimesteps, args["timewidth"])
chanblocks = cld(mset.nchans, args["chanwidth"])
@info "Calibrating with $timeblocks independent timeblocks and $chanblocks independent chanblocks"

# Initialize final Jones matrix to Identity
# [pol, pol, antid, channels, timesteps]
jones = zeros(Complex{Float64}, 4, mset.nants, chanblocks, timeblocks)
jones[1, :, :, :] .= 1
jones[4, :, :, :] .= 1
converged = zeros(Bool, chanblocks, timeblocks)
@debug "Created Jones array: $(size(jones))"

# Initialise channel for sending work from the producer to consumers.
# The producer is responsible for loading data from disk (or predicting it) and
# distributing the chunks via the channel to the consumers, who do the actual calibration.
# Since loading data is IO bound, we only have one producer, but we inialise as many
# consumers as we have threads, since this is CPU bound.
ch = Channel{
    Tuple{Array{ComplexF32, 3}, Array{ComplexF32, 3}, Array{Float32, 3}, Vector{Int}, Vector{Int}, Int, Int}
}(args["nbatch"])

# Initialise consumer threads
tasks = Task[]
for _ in 1:Threads.nthreads()
    task = Threads.@spawn consumer(ch, jones, converged, args)
    push!(tasks, task)
end

# And run producer in this thread
producer(ch, mset, comps, beam, args)

# Wait on consumers to complete their calibration
close(ch)  # When the channel is drained, this signals to workers to return
for task in tasks
    wait(task)
end

# Invert jones Matrices so that they can be applied as J D J^H
for timeblock in axes(jones, 4), chanblock in axes(jones, 3), antid in axes(jones, 2)
    @views invA!(jones[:, antid, chanblock, timeblock])
end

save(File(format"JLD", args["solution"]), "jones", jones)

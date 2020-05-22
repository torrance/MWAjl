using ArgParse
using Logging
using MWAjl

function main(args)
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
    if args["apply-beam"] && mset.mwadelays !== nothing
        global beam = Beam(mset.mwadelays)
    elseif mset.mwadelays == nothing
        @error "Unable to load MWA tile delays from measurement set, which is required with --apply-beam"
        exit(1)
    else
        global beam = nothing
    end

    # Set channel width if --chanwidth==0
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

    # Set default --nbatch if --nbatch==0
    if args["nbatch"] < 1
        args["nbatch"] = max(1, floor(Int, 200 / args["chanwidth"]))
        @info "Set --nbatch = $(args["nbatch"])"
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

    # Write out the solution
    open(args["solution"], "w") do f
        # First expand the Jones array for each channel
        expandedjones = zeros(Complex{Float64}, 4, mset.nants, mset.nchans, timeblocks)
        for chan in 1:mset.nchans
            chanblock = ceil(Int, chan / args["chanwidth"])
            expandedjones[:, :, chan, :] .= jones[:, :, chanblock, :]
        end
        # Now write the solution to the binary format we've inherited from AOCal
        writesolution(f, expandedjones, mset.unique_timesteps[1], mset.unique_timesteps[end])
    end

    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    s = ArgParseSettings("Calibrate Murchison Widefield Array (MWA) observations. See https://torrance.github.io/MWAjl/")
    @add_arg_table! s begin
        "--model", "-m"
            help="The path to the model file (aoskymodel 1.2 format). If absent, will use the --modelcolumn column (e.g. during self-calibration)."
            arg_type=String
        "--chanwidth", "-c"
            help="Find calibration solutions for blocks of channels of --chanwidth. Set to 0 to find a single solution for all channels. Larger values may provide more signal and aid in finding a better calibration solution but at the expense of ignoring frequency-dependenct changes to the calibration solution."
            arg_type=Int
            default=1
        "--timewidth", "-t"
            help="Find calibrations solutions for blocks of timesteps of --timewidth. Defaults to 0, which implies an infinite time width. The duration of a timestep depends on the resolution of the Measurement Set and any averaging in time that may have be performed. e.g. A 2 minute observation with time resolution 4 s will have 30 timesteps; setting --timewidth 10 will result in 3 independent solutions in time."
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
            help="Determines whether a particular calibration solution has converged. Consists of two values: mininum accuracy and stopping accuracy. These two values determine whether 1) a solution has sufficiently converged that we can accept its answer and 2) whether it has converged enough that we may stop prior to reaching --max-iterations. These values are tested after each iteration by comparing the magnitude difference between the current solution and prior solution."
            arg_type=Float64
            nargs=2
            default=[1E-5, 1E-8]
        "--max-iterations", "-i"
            help="The maximum number of iterations allowed when solving for a single solution unit (i.e. for a given channel and time block). Usually a good solution can be found in 10-20 iterations. If you consistently hit the default limit, consider relaxing the --tolerance stopping accuracy. Higher --max-iteration values than the default are of dubious benefit."
            arg_type=Int
            default=50
        "--apply-beam"
            help="Apply the MWA beam during model prediction to correct model flux values from their true values to their apparent values. The beam is calculated in full polarization for each timestep and for each coarse channel."
            action=:store_true
        "--datacolumn"
            help="The uncalibrated data column in the Measurement Set. Case sensitive."
            arg_type=String
            default="DATA"
        "--modelcolumn"
            help="The model data column in the Measurement Set that is used when --model is not provided. Case sensitive."
            arg_type=String
            default="MODEL_DATA"
        "--nbatch"
            help="When reading from --datacolumn and --modelcolumn, this parameter specifies how many --chanwidth by --timewidth arrays to read. This parameter only affects performance and memory usage. Ideally, this needs to be large enough so that calibration workers aren't waiting on data to be read from disk, but not so large that too long is spent waiting at the start of the program for the calibration workers to begin work. Reducing this value will reduce maximum memory usage proportionally. Default (0) sets this to 200 / --chanwidth."
            arg_type=Int
            default=0
        "--gpu"
            help="Use GPU acceleration for model prediction. Requires a CUDA capable graphics card and associated CUDA drivers."
            action=:store_true
        "--verbose"
            help="Print logging information about what calibrate is doing."
            action=:store_true
        "--debug"
            help="Print lots more logging information about what calibrate is doing."
            action=:store_true
        "mset"
            help="The path to the CASA Measurement Set to be calibrated."
            required=true
        "solution"
            help="The output solution file name and path. e.g. folder/to/go/solutions.bin"
            required=true
    end

    args = parse_args(s)
    main(args)
end
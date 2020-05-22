# MWAjl

A Julia library for calibrating measurement sets produced by the Murchison Widefield Array (MWA).

This is a rewrite of the existing calibration software using a near-identical algorithm (known variously as 'Mitchcal' or 'Stefcal'), but optimised for speed and to take advantage of graphics processing units (GPUs) if available.

```@contents
Pages = ["index.md"]
Depth = 5
```

## Installation

The package has been tested in Julia 1.4 and can be installed from the Julia package manager. Start the Julia REPL and press `]` to enter the package manager, then run `add https://github.com/torrance/MWAjl`, eg:

    $ julia
    julia> ]
    pkg> add https://github.com/torrance/MWAjl

There are additional build steps required to link against Casacore and the MWA beam.

It is recommened to install PackageCompiler to ensure fast start-up times when using `calibrate.jl`.

## Usage

### Overview

There are two modes of running calibration: from a sky model or self-calibration. These two modes differ only in how they obtain the model visibilities. In the sky model mode, the visibilities are predicted from the `--model` file, whilst the self-calibration mode uses the existing `MODEL_DATA` column as its model (which has usually been populated during an earlier imaging round).

Self-calibration is the default mode, unless a sky model is provided by passing the `--model` parameter.

### Options

```
usage: calibrate.jl [-m MODEL] [-c CHANWIDTH] [-t TIMEWIDTH]
                    [--minuv MINUV] [--maxuv MAXUV]
                    [-a TOLERANCE TOLERANCE] [-i MAX-ITERATIONS]
                    [--apply-beam] [--datacolumn DATACOLUMN]
                    [--modelcolumn MODELCOLUMN] [--nbatch NBATCH]
                    [--gpu] [--verbose] [--debug] [-h] mset solution

Calibrate Murchison Widefield Array (MWA) observations. See
https://torrance.github.io/MWAjl/

positional arguments:
  mset                  The path to the CASA Measurement Set to be
                        calibrated.
  solution              The output solution file name and path. e.g.
                        folder/to/go/solutions.bin

optional arguments:
  -m, --model MODEL     The path to the model file (aoskymodel 1.2
                        format). If absent, will use the --modelcolumn
                        column (e.g. during self-calibration).
  -c, --chanwidth CHANWIDTH
                        Find calibration solutions for blocks of
                        channels of --chanwidth. Set to 0 to find a
                        single solution for all channels. Larger
                        values may provide more signal and aid in
                        finding a better calibration solution but at
                        the expense of ignoring frequency-dependenct
                        changes to the calibration solution. (type:
                        Int64, default: 1)
  -t, --timewidth TIMEWIDTH
                        Find calibrations solutions for blocks of
                        timesteps of --timewidth. Defaults to 0, which
                        implies an infinite time width. The duration
                        of a timestep depends on the resolution of the
                        Measurement Set and any averaging in time that
                        may have be performed. e.g. A 2 minute
                        observation with time resolution 4 s will have
                        30 timesteps; setting --timewidth 10 will
                        result in 3 independent solutions in time.
                        (type: Int64, default: 0)
  --minuv MINUV         Exclude baselines shorter than this length
                        from affecting calibration solution (metres).
                        (type: Float64, default: 0.0)
  --maxuv MAXUV         Exclude baselines longer than this length from
                        affecting calibration solutions (metres).
                        (type: Float64, default: 9.0e99)
  -a, --tolerance TOLERANCE TOLERANCE
                        Determines whether a particular calibration
                        solution has converged. Consists of two
                        values: mininum accuracy and stopping
                        accuracy. These two values determine whether
                        1) a solution has sufficiently converged that
                        we can accept its answer and 2) whether it has
                        converged enough that we may stop prior to
                        reaching --max-iterations. These values are
                        tested after each iteration by comparing the
                        magnitude difference between the current
                        solution and prior solution. (type: Float64,
                        default: [1.0e-5, 1.0e-8])
  -i, --max-iterations MAX-ITERATIONS
                        The maximum number of iterations allowed when
                        solving for a single solution unit (i.e. for a
                        given channel and time block). Usually a good
                        solution can be found in 10-20 iterations. If
                        you consistently hit the default limit,
                        consider relaxing the --tolerance stopping
                        accuracy. Higher --max-iteration values than
                        the default are of dubious benefit. (type:
                        Int64, default: 50)
  --apply-beam          Apply the MWA beam during model prediction to
                        correct model flux values from their true
                        values to their apparent values. The beam is
                        calculated in full polarization for each
                        timestep and for each coarse channel.
  --datacolumn DATACOLUMN
                        The uncalibrated data column in the
                        Measurement Set. Case sensitive. (default:
                        "DATA")
  --modelcolumn MODELCOLUMN
                        The model data column in the Measurement Set
                        that is used when --model is not provided.
                        Case sensitive. (default: "MODEL_DATA")
  --nbatch NBATCH       When reading from --datacolumn and
                        --modelcolumn, this parameter specifies how
                        many --chanwidth by --timewidth arrays to
                        read. This parameter only affects performance
                        and memory usage. Ideally, this needs to be
                        large enough so that calibration workers
                        aren't waiting on data to be read from disk,
                        but not so large that too long is spent
                        waiting at the start of the program for the
                        calibration workers to begin work. Reducing
                        this value will reduce maximum memory usage
                        proportionally. Default (0) sets this to 200 /
                        --chanwidth. (type: Int64, default: 0)
  --gpu                 Use GPU acceleration for model prediction.
                        Requires a CUDA capable graphics card and
                        associated CUDA drivers.
  --verbose             Print logging information about what calibrate
                        is doing.
  --debug               Print lots more logging information about what
                        calibrate is doing.
  -h, --help            show this help message and exit
```

### When to use `--minuv` and `--maxuv`

These two parameters allow you to exclude baselines by their length (metres) from being used to form calibration solutions. That length is the distance between the antennas in 3D space; it is not, for example, a projection onto a _uv_ plane based on the current phase centre, and so will not take into account foreshortening effects.

The primary reason to use these parameters is when using an incomplete model of the sky (which is always the case!) that might be insensitive to large scale emission features or might be lacking in resolution.

For example the GLEAM catalogue is commonly used to calibrate MWA observations. The catalogue does not include large scale Galactic emission or supernova remnants. Since these features are of large angular extent, most of their power will be detected on the short baselines. On the other hand, GLEAM used the MWA Phase I array which had a maximum baseline length of just a few kilometres; some of the apparent 'point sources' in GLEAM may actually be resolved to multiple sources in the Phase II array. We can therefore set `--minuv 100` and `--maxuv 2600` (as an example only!) to use baselines in this interval since it is only in this range that our skymodel is fairly accurate.

On the other hand, these parameters will exclude data and so make it harder for `calibrate` to find a good solution amongst all the visibility noise.

Note also that one set of values might be appropriate when calibrating using a skymodel (e.g. GLEAM), but that another set of values might be appropriate when self-calibrating. In the latter case, it would usually make sense to not impose a `--maxuv` cutoff at all.

### Understanding `--tolerance`

The `--tolerance` parameter is used to define the convergence of a solution. It is a measure of how much the solution is changing between calibration iterations, and it controls when `calibrate` thinks what it's found has settled down to a stable value and is 'good enough'. A smaller tolerance indicates a higher degree of required stability of the solution.

The parameter accepts two values: `min-accuracy` and `stopping-accuracy`. If `calibrate` reaches a stablity between calibration iterations that is smaller than `stopping-accuracy`, it stops immediately and accepts the solution. However, if by the time `calibrate` reaches the `--max-iterations` parameter and the stability is less than `min-accuracy`, then we will also accept this as a valid solution. If the solution hasn't settled down even to this level of stability, we mark it as failed.

It's important to note there is a tradeoff between the performance of `calibrate` and `--tolerance`. Very low values of `--tolerance` may force `calibrate` to consistently hit `--max-iterations`, without appreciably increasing the quality of the calibration solution.

For initial calibration (i.e. using `--model`, which is only a rough and incomplete skymodel), relaxed values of tolerance such as `--tolerance 1E-5 1E-6` are recommended. Self-calibration typically can reach much higher stability quite quickly, and values such as `--tolerance 1E-5 1E-8` are good.

If you're consistently hitting the `--max-iterations` threshold, consider whether your `--tolerance` is too low. A good solution can usually be reached in only 10-20 iterations.

The stability value inherits its definition from the original `calibrate` and is calculated as follows. Each iteration of the Mitchcal/Stefcal algorithm provides a new calibration solution ``J[4, \text{NANTS}]``. This solution is a 4 x NANTS matrix, where we have a solution for each of the NANTS (e.g. 128) antennas, and for each of the 4 instrumental polarisations. Now, we take the prior solution ``J_1`` and next proposed solution ``J_2`` and we calculate the squared Eucliden distance (elementwise) between the two solutions, ``\Delta = ||J_1 - J_2||^2``. Then we find the mean value of this distance for each of the four instrumental polarisations across all antennas (i.e. the mean of each row of ``\Delta``). Finally, the stability is defined as the maximum of these four values. That is the stability ``S`` is defined as (using 1-indexed notation):

```math
S := \text{max}(\text{mean}(\Delta, \text{axis}=2))
```

## Managing memory usage

There are three parameters that affect how `calibrate` uses memory. When `calibrate` runs, it independently calibrates chunks of visibility data that are `--timewidth` long in duration (measured in timesteps, not seconds) and `--chanwidth` wide in frequency space (measured in channels). And `calibrate` reads from disk `--nbatch` chunks of this data at once.

`--nbatch` is the first option you should change if you want to reduce memory usage. It defaults to `--nbatch 200` which means it reads in 200 chunks of visibility data that are `--timewidth` by `--chanwidth` in size. You can reduce memory usage by half by simply setting `--nbatch 100`. Note that `--nbatch` is a performance optimization: it speeds up performance to read more than just a single channel at a time; it is also faster to send a large array of data to the GPU for prediction rather than multiple small arrays; and prefetching lots of data helps keep the calibration workers from becoming idle. If `--nbatch` gets so small that calibration workers routinely end up idle, then performance will start to suffer.

The other two parameters are `--timewidth` and `--chanwidth`. The default values are `--timewidth 0` which corresponds to all timesteps and `--chanwidth 1` which corresponds to a single channel. However, it may be desirable to increase signal by setting `--chanwidth 2`, for example, and this will double the size of the data in memory. Or, perhaps, you think that a calibration solution varies across the 30 timesteps and instead choose to set `--timewidth 10` so that you have 3 separate solutions in time; this will reduce memory usage by 3 times. If you're feeling crazy, setting `--timewidth 0` and `--chanwidth 0` will force `calibrate` to attempt to find a single solution for all antennas in time and frequency, and will result in the entire Measurement Set being read into memory.


## Performance notes

### Expected runtime

On my own (slightly underpowered, 4/8 core, 64 Gb memory) desktop, the following table gives a rough idea of how long it _should_ take to calibrate a 2 minutes MWA observation (128 tiles):

| GPU | Disk type | Components | Time   |
| :-- | :-------- | :--------- | :----- |
| Yes | SSD       | 201        | 1m40s  |
| Yes | SSD       | 1028       | 4m40s  |
| No  | SSD       | 201        | 6m36s  |
| No  | SSD       | 1028       | 24m50s |
| n/a | SSD       | selfcal    | 1m40s  |
| n/a | HDD       | selfcal    | 2m30s  |

Of course, these are only approximate values for my own desktop machine and are a snapshot in time. But they should indicate what you should roughly expect, and whether you're encountering slower than expected results.

### Model sources

As can be seen above, performance suffers as the number of sources in `--model` increases. This is because each source (actually, each component - a single source may be comprised of multiple components) has to be added to the visibilities by way of direct Fourier transform. The computational cost of visibility prediction scales linearly with the number of sources.

Often, it is not necessary to have a large number of sources. The vast majority of flux in the MWA field of view is often produced by just a handful of sources, and it's usually quite sufficient to calibrate with a model of just 100-200 sources to arrive at a good, initial calibration. The second self-calibration step can take you the rest of the way.

### GPU acceleration

GPU acceleration is used to do the prediction step, that is, turning a `--model` file into a set of model visibilities. The performance benefit of this varies between computers, but on my own computer I see ~6x speed up by using the GPU as opposed to my 4 cores (8 hyperthreads). This factor also depends on how many components make up your model.

Note that for self-calibration the GPU is not used at all. As it turned out, the actual calibration step is currently very fast on a CPU, very awkward to perform on a GPU due to a reduction step, and by far the least slowest part of running calibrate.

### IO issues

A significant bottleneck of `calibrate` is reading data from disk. `calibrate` uses the CasaCore library to read visibility data that is in the Measurement Set format. The `--nbatch` option controls how we read data from the Measurement Set, but usually we read a handful of channels at a time (as opposed all channels) to reduce overall memory usage, and this data is not contiguous on disk. For a hard disk drive (i.e. spinning platters), the extra time associated with seeking can add overhead. And on non-local storage such as in a supercomputing cluster, the time spent reading from disk may be the most significant factor in the performance of `calibrate`.

Additionally, even when the Measurement Set is stored entirely in memory, there is still some overhead associated with CasaCore itself. For example, on my system it takes 5 s to read the first 100 channels from a Measurement Set that is stored in a tmpfs. It is certainly worth investigating with alternative storages such as HDF5 suffer from the same overhead.

## Todo:

* explain --tolerance, and recommended settings
* memory requirements, and mitigations with --timewidth and --chanwidth and --nbatch
* what does GPU acceleration do?
* document stefcal
* format of solution.bin files
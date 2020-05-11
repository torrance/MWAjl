# MWAjl

A Julia library for calibrating measurement sets produced by the Murchison Widefield Array (MWA).

This is a rewrite of the existing calibration software using a near-identical algorithm (known variously as 'Mitchcal' or 'Stefcal'), but optimised for speed and to take advantage of graphics processing units (GPUs) if available.

## Installation

The package has been tested in Julia 1.3 and can be installed from the Julia package manager. Start the Julia REPL and press `]` to enter the package manager, then run `add https://github.com/torrance/MWAjl`, eg:

    $ julia
    julia> ]
    pkg> add https://github.com/torrance/MWAjl

There are additional build steps required to link against Casacore and the MWA beam.

It is recommened to install PackageCompiler to ensure fast start-up times when using `calibrate.jl`.

## Usage

There are two modes of running calibration: from a sky model or self-calibration. These two modes differ only in how they obtain the model visibilities. In the sky model mode, the visibilities are predicted from the `--model` file, whilst the self-calibration mode uses the existing `MODEL_DATA` column as its model (which has usually been populated during an earlier imaging round).

Self-calibration is the default mode, unless a sky model is provided by passing the `--model` parameter.
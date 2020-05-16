struct MeasurementSet
    table::Table
    nrows::Int
    nants::Int
    nchans::Int
    ntimesteps::Int
    unique_timesteps::Array{Float64, 1}
    freqs::Array{Float32, 1}
    phasecenter::Position
    mwadelays::Union{Array{Int32, 1}, Nothing}
end

function MeasurementSet(path::String)
    mset = Table(path)

    # Timesteps
    time = column(mset, "TIME")
    unique_timesteps = sort(unique(time))
    ntimesteps = length(unique_timesteps)

    # Rows
    nrows = length(time)

    # Frequency
    spw = Table(path * "/SPECTRAL_WINDOW")
    freqs = column(spw, "CHAN_FREQ")
    if size(freqs, 2) > 1
        throw(ErrorException("There is more than one spectral window present in the measurement set. We only know how to handle one."))
    end
    freqs = freqs[:, 1]
    nchans = length(freqs)

    # Antennas
    antennas = Table(path * "/ANTENNA")
    nants = size(column(antennas, "POSITION"), 2)

    # Phase center
    pos0 = Position(
        column(Table(path * "/FIELD"), "PHASE_DIR")[:, 1]...
    )

    # MWA tile delays
    mwadelays = column(Table(path * "/MWA_TILE_POINTING"), "DELAYS")[:, 1]

    @info """
        Measurement set: $path
            Rows: $nrows
            Antennas: $nants
            Timesteps: $ntimesteps
            Channels: $nchans ($(freqs[1] / 1E6) - $(freqs[end] / 1E6) MHz)
            Phase center RA = $(rad2deg(pos0.ra)) Dec = $(rad2deg(pos0.dec))
            MWA delays: $mwadelays"""

    return MeasurementSet(
        mset, nrows, nants, nchans, ntimesteps, unique_timesteps, freqs, pos0, mwadelays
    )
end

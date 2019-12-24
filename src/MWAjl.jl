module MWAjl

include("skymodel.jl")
include("skymodel_parser.jl")
export parse_model, Source, Component, SED, Measurements

include("casacore.jl")
export Table, CasaError, column

include("stefcal.jl")
export calibrate

end # module

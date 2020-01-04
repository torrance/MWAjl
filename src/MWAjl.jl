module MWAjl

include("beam.jl")
export Beam, beamjones

include("skymodel.jl")
include("skymodel_parser.jl")
export parse_model, Source, Component, SED, Measurements, Position, lmn

include("casacore.jl")
export Table, CasaError, column, taql, Frame, radec2altaz

include("predict.jl")
export predict

include("stefcal.jl")
export calibrate!, innerloop

include("matrix2x2.jl")
export Matrix2x2toArray4, Array4toMatrix2x2, AxB!, AxBH!, AHxB!, plusAxB!, plusAxBH!, plusAHxB!, AdivB!, invA!

include("utils.jl")
export sanitize!

const mwapb = PyNULL()
function __init__()
    # Python modules must be loaded as part of the module init, otherwise
    # empty NULL pointers will be saved following module precompilation.
    copy!(mwapb, pyimport("mwa_pb.primary_beam"))
end

end

module MWAjl

include("beam.jl")
export Beam, calc_jones, closest_freq

include("skymodel.jl")
include("skymodel_parser.jl")
export parse_model, gettoken, stokes, Source, Component, SED, Measurements, Position, lmn, instrumental

include("components/gaussian.jl")
export Gaussian, predictgaussian!

include("casacore.jl")
export Table, CasaError, column, taql, Frame, radec2altaz

include("measurementset.jl")
export MeasurementSet

include("predict.jl")
export predict

include("solutionfile.jl")
export writesolution

include("stefcal.jl")
export calibrate!, calibrationloop

include("matrix2x2.jl")
export Matrix2x2toArray4, Array4toMatrix2x2, AxB!, AxBH!, AHxB!, plusAxB!, plusAxBH!, plusAHxB!, AdivB!, invA!

include("workers.jl")
export consumer, producer

include("utils.jl")
export sanitize!

end

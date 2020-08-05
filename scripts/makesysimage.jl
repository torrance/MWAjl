using Pkg
using Pkg.Artifacts
using PackageCompiler

Pkg.activate("..")

project = joinpath(@__DIR__, "..")
calibrate = joinpath(@__DIR__, "../bin/calibrate.jl")
testdata = artifact"testdata"

temp_dir = mktempdir()
cpu_precompile = joinpath(temp_dir, "cpu-precompile.jl")
gpu_precompile = joinpath(temp_dir, "gpu-precompile.jl")
precompile = joinpath(temp_dir, "precompile.jl")


run(`julia --trace-compile=$cpu_precompile --project=$project $calibrate --model $testdata/data/model.txt --apply-beam -a 1E-5 1E-6 --verbose $testdata/data/1217094992-ch8-19.ms /dev/null`)

# Currently, CUDA drivers are not available during Singulairty build processes.
# run(`julia --trace-compile=$gpu_precompile --project=$project $calibrate --model $testdata/data/model.txt --apply-beam -a 1E-5 1E-6 --gpu --verbose $testdata/data/1217094992-ch8-19.ms /dev/null`)
run(`touch $gpu_precompile`)

open(precompile, "w") do io
    write(
        io,
        read(`sort -u $cpu_precompile $gpu_precompile`, String),
    )
end

# Remove MWAjl as for some reason this causes the sysimage to segfault. This needs to be fixed.
PackageCompiler.create_sysimage([
    :ArgParse, :CUDAdrv, :CUDAnative, :CuArrays, :HDF5, :LinearAlgebra, :Pkg, :StaticArrays, :Statistics
], sysimage_path=joinpath(project, "libs/MWAjl.so"), precompile_statements_file=precompile)

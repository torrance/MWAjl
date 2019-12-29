using MWAjl
using BenchmarkTools

nants = 128
ntimesteps = 30
nchans = 4
nrows = Int((nants * (nants - 1)) / 2 * ntimesteps)
println("Rows: $nrows")

data = zeros(ComplexF32, 2, 2, nchans, nrows)
model = zeros(ComplexF32, 2, 2, nchans, nrows)
ants1 = zeros(Int, nrows)
ants2 = zeros(Int, nrows)
weights = randn(Float32, size(data)...)
weights .= abs.(1 .+ weights)

# Initialize data
data[1, 1, :, :] .= 1
data[2, 2, :, :] .= 1
model[1, 1, :, :] .= 1
model[2, 2, :, :] .= 1

# Create true jones
jones = (reshape([1 0; 0 1], 2, 2, 1) .+ 0.2 * randn(Float32, 2, 2, 128)) .* exp.(2im * π * rand(Float32, 2, 2, 128))

# Uncalibrate data
row = 1
for timestep in 1:ntimesteps, ant1 in 1:nants, ant2 in 1:nants
    global row
    if ant1 >= ant2
        continue
    end
    ants1[row] = ant1
    ants2[row] = ant2
    for chan in 1:nchans
        data[:, :, chan, row] = jones[:, :, ant1] * data[:, :, chan, row] * jones[:, :, ant2]'
    end
    row += 1
end

# Find solution
jones1 = similar(jones)
jones1 .= 0
jones1[1, 1, :] .= 1
jones1[2, 2, :] .= 1

top = similar(jones1)
bot = similar(jones1)
z = zeros(ComplexF64, 2, 2)

t = @benchmark innerloop($data, $model, j, $ants1, $ants2, $top, $bot, $z) setup=(j = copy($jones1)) evals=1 samples=50 seconds=999
show(stdout, MIME"text/plain"(), t)

t = @benchmark calibrate!(j, $data, $model, $weights, $ants1, $ants2, 50, 1E-5, 1E-8) setup=(j = copy($jones1)) evals=1 samples=5 seconds=999
show(stdout, MIME"text/plain"(), t)

# @code_warntype calibrate!(zeros(ComplexF64, 2, 2, 128), zeros(ComplexF32, 2, 2, 4, 1000), zeros(ComplexF32, 2, 2, 4, 1000), zeros(Float32, 2, 2, 4, 1000), zeros(Int, 1000), zeros(Int, 1000), 50, 1E-5, 1e-8)
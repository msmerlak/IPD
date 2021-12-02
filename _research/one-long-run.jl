

using Distributed
using DataFrames, CSV#, Random, CSV, Dates

#PARALLEL = true
PROCESSES = 10#Sys.CPU_THREADS - 1
addprocs(PROCESSES)

@everywhere begin
    import Pkg
    Pkg.activate(".")
end
@everywhere using DrWatson

mkdir(plotsdir("mixed", string(today())))
mkdir(datadir("mixed", string(today())))


@everywhere begin
    @quickactivate "IPD"
    include("../src/utils.jl")
    include("../src/IPD.jl")
    include("../src/constants.jl")
end

include(srcdir("plotting.jl"))

### distributions

σ, n, t = 5e-3, 500, 10
init = rand(4)
multiplicative = false
T = 1_000_000

p = Dict(:RSTP => RSTP, :n => n, :t => t, :σ => σ, :multiplicative => multiplicative, :ϵ => 1e-10, :init_share => false, :m => 0.)

# model = create_model(p, rng = MersenneTwister())
models = [create_model(p; rng = MersenneTwister(rand(UInt)), initial_strategy = rand(4)) for _ in 1:PROCESSES] 

adata, mdata = ensemblerun!(models, player_step!, WF_sampling!, T, adata = [(:fitness, mean)], parallel = true)

CSV.write(datadir("mixed", "one-long-run.csv"), adata)


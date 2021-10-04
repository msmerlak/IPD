using DrWatson

using Distributed
using DataFrames, MultivariateStats, StatsBase, CSV
include(srcdir("plotting.jl"))

PARALLEL = true
PROCESSES = Sys.CPU_THREADS-1
PARALLEL && addprocs(PROCESSES)
@everywhere begin
    include("../src/IPD.jl")
end



models = [create_model(Dict(:RSTP => ID_PAYOFFS, :n => 200, :t => 10, :σ => 5*1e-3, :multiplicative => true, :init => rand(4), :seed => seed, :ϵ => 1e-10)) for seed in rand(UInt8, PROCESSES)]

adata, _ = ensemblerun!(models, player_step!, WF_sampling!, 10000, adata = [(:score, mean), (:strategy, mean)], parallel = true)

CSV.write(datadir("mixed", savename(models[1].properties, "csv")), adata)

plot_probs(adata)

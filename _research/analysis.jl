using DrWatson

using Distributed
using DataFrames, MultivariateStats, StatsBase, CSV
include(srcdir("plotting.jl"))
#include(srcdir("ensemblerun.jl"))

PARALLEL = true
PROCESSES = 15
PARALLEL && addprocs(PROCESSES)
@everywhere begin
    include("../src/IPD.jl")
end

p = Dict(:RSTP => ID_PAYOFFS, :n => 500, :t => 10, :σ => 5*1e-2, :multiplicative => false, :ϵ => 1e-10, :init_share => false, :m => 0.)

models = [create_model(p; rng = MersenneTwister(rand(UInt)), initial_strategy = rand(4)) for _ in 1:PROCESSES] 

adata, _ = ensemblerun!(models, player_step!, WF_sampling!, 5000, adata = [(:fitness, mean), (:strategy, mean)], parallel = true)

CSV.write(datadir("mixed", savename(models[1].properties, "csv")), adata)

plot_probs(adata)

begin
    p = plot()
    for ensemble in groupby(adata, :ensemble)
        plot!(p, ensemble.mean_share)
    end
    current()
end

using Distributed
using DataFrames, CSV, Dates, Random, RandomNumbers

#PARALLEL = true
PROCESSES = Sys.CPU_THREADS - 1
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
    include("../src/IPD.jl")
    include("../src/utils.jl")
    include("../src/constants.jl")
end

#include(srcdir("plotting.jl"))
using StatsPlots, StatsBase, Query, Statistics

@everywhere function initialize(; n, σ)
    p = Dict(:RSTP => RSTP, :n => n, :σ => σ)
    return create_model(p; initial_strategy = rand(4))
end


### varying N

params = Dict(:n => collect(20:20:100), :σ => 1e-2)

adata, _ = paramscan(params, initialize; 
    agent_step! = player_step!, 
    model_step! = WF_sampling!, 
    n = 1_000_000, 
    adata = [(:fitness, mean), (:strategy, mean)], 
    parallel = true
    )


@df adata plot(
    :step,
    :mean_fitness,
    group = :n
    )

plts = []
for n in params[:n]
    push!(plts, 
        plot(adata[adata.n .== n, :mean_fitness],
        ylims = (.8, 3), label = n)
    )
end

plot(plts..., layout = (5,1), dpi= 500)
savefig(plotsdir("vary-N.png"))
    
movingaverage(g, n = 1000) = [i < n ? mean(g[begin:i]) : mean(g[i-n+1:i]) for i in 1:length(g)]

df = adata |> @filter(_.n == 80) |> DataFrame
plot(
    movingaverage(df.mean_fitness, 10000)
)


function sojourn_times(f)
    smoother = movingaverage(f, 10000)
    x, y = rle(smoother .> 2)
    return y[x]
end


plot(
    map(n -> mean(sojourn_times( adata[adata.n .== n, :mean_fitness] )), params[:n])
)

### varying σ

params = Dict(:n => 50, :σ => collect(.005:0.01:0.05))

adata_mut, _ = paramscan(params, initialize; 
    agent_step! = player_step!, 
    model_step! = WF_sampling!, 
    n = 1_000_000, 
    adata = [(:fitness, mean), (:strategy, mean)], 
    parallel = true
    )

plts = []
for σ in params[:σ]
    push!(plts, 
        plot(adata_mut[adata_mut.σ .== σ, :mean_fitness],
        ylims = (.8, 3), label = σ)
    )
end
plot(plts..., layout = (length(params[:σ]),1), dpi= 300)
savefig(plotsdir("vary-σ.png"))

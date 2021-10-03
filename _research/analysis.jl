using DrWatson
@quickactivate 
using Distributed
using DataFrames, MultivariateStats, StatsBase
include(srcdir("plotting.jl"))

PARALLEL = true

PARALLEL && addprocs(Sys.CPU_THREADS)
@everywhere include("src/IPD.jl")

models = [create_model(
        mutational_effect = 5*1e-3, 
        tournament_size = 10, 
        popsize = 200, 
        multiplicative = false, 
        seed = x,
        initial_strategy = TFT
        ) for x in rand(UInt8, 20)]

adata, _ = ensemblerun!(models, player_step!, WF_sampling!, 5000, adata = [(:score, mean), (:strategy, mean)], parallel = true)

# plot([mean(df.score) for df in groupby(filter(player -> length(player.score)> 0, adata), :step)], xlabel = "generation", ylabel = "mean score")

plot_probs(adata)

# p = plot()
# for realization in groupby(adata, :ensemble)
#     st = [mean(df.strategy) for df in groupby(realization, :step)]
#     mean_strategies = reduce(hcat, st) |> transpose 

#     plot!(p, mean_strategies[:, 1], mean_strategies[:, 2], mean_strategies[:, 3])
# end
# current()

p = plot()
for realization in groupby(adata, :ensemble)
    plot!(reduce(hcat, realization.mean_strategy)[4, :], legend = false)
end
current()



# using ImageCore, CairoMakie
# f = colorsigned() ∘ scalesigned(0., 0, 100.0)
# abm_plot(model, ac = a -> f(a.score))[1]

## not parallel

# model = create_model(multiplicative = true, seed = rand(UInt8))
# adata, _ = run!(model, player_step!, WF_sampling!, 1000, adata = [:score, :strategy])


# result =  groupby(adata, :step)
# st = [mean(df.strategy) for df in result]
# mean_strategies = reduce(hcat, st) |> transpose 

# plot(mean_strategies)


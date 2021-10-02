PARALLEL = true

using Distributed
using DataFrames, Plots, InteractiveDynamics, MultivariateStats, StatsBase
include("../src/plotting.jl")

PARALLEL && addprocs(10)
@everywhere include("../src/IPD.jl")


models = [create_model(
        mutational_effect = 1e-2, 
        tournament_size = 10, 
        popsize = 200, 
        multiplicative = false, 
        seed = x,
        initial_strategy = TFT
        ) for x in rand(UInt8, 10)]

adata, _ = ensemblerun!(models, player_step!, WF_sampling!, 1000, adata = [:score, :strategy], parallel = true)


plot([mean(df.score) for df in groupby(adata, :step)], xlabel = "generation", ylabel = "mean score")


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
    st = [mean(df.strategy) for df in groupby(realization, :step)]
    mean_strategies = reduce(hcat, st) |> transpose 
    #println(pca)
    plot!(mean_strategies[:, 4])
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

using DrWatson
@quickactivate "IPD"

using Distributed
addprocs(10)
@everywhere include("IPD.jl")

models = [create_model(
        mutational_effect = 1e-2, 
        tournament_size = 5, 
        popsize = 500, 
        multiplicative = false, 
        ϵ = 1e-50, seed = x,
        initial_strategy = TFT
        ) for x in rand(UInt8, 5)]

adata, _ = ensemblerun!(models, player_step!, WF_sampling!, 500, adata = [:score, :strategy], parallel = true)

using DataFrames, Plots, InteractiveDynamics

result =  groupby(adata, :step)



using MultivariateStats
p = plot()
for realization in groupby(adata, :ensemble)
    st = [mean(df.strategy) for df in groupby(realization, :step)]
    mean_strategies = reduce(hcat, st) |> transpose 
    pca = projection(fit(PCA, mean_strategies, maxoutdim = 2))
    #println(pca)
    plot!(p, pca[:, 1], pca[:, 2])
end
current()




# mm = model.multiplicative
# plot(mean_strategies, label = ["CC" "CD" "DC" "DD"], title = "multiplicative = $mm")


# using ImageCore, CairoMakie
# f = colorsigned() ∘ scalesigned(0., 0, 100.0)
# abm_plot(model, ac = a -> f(a.score))[1]
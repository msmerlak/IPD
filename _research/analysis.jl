PARALLEL = true

using Distributed
PARALLEL && addprocs(10)
@everywhere include("IPD.jl")

models = [create_model(
        mutational_effect = 5*1e-3, 
        tournament_size = 10, 
        popsize = 100, 
        multiplicative = false, 
        seed = x,
        initial_strategy = rand(4)
        ) for x in rand(UInt8, 10)]

adata, _ = ensemblerun!(models, player_step!, sampling!, 1000, adata = [:score, :strategy], parallel = true)


using DataFrames, Plots, InteractiveDynamics
result =  groupby(adata, :step)


using MultivariateStats, StatsBase
p = plot()
for realization in groupby(adata, :ensemble)
    st = [mean(df.strategy) for df in groupby(realization, :step)]
    mean_strategies = reduce(hcat, st) |> transpose 
    #pca = projection(fit(PCA, mean_strategies, maxoutdim = 2))
    #println(pca)
    plot!(p, mean_strategies[:, 1], mean_strategies[:, 2])
end
current()


p = plot()
for realization in groupby(adata, :ensemble)
    st = [mean(df.strategy) for df in groupby(realization, :step)]
    mean_strategies = reduce(hcat, st) |> transpose 
    #pca = projection(fit(PCA, mean_strategies, maxoutdim = 2))
    #println(pca)
    plot!(p, mean_strategies[:, 3], mean_strategies[:, 4])
end
current()

p = plot()
for realization in groupby(adata, :ensemble)
    st = [mean(df.strategy) for df in groupby(realization, :step)]
    mean_strategies = reduce(hcat, st) |> transpose 
    #pca = projection(fit(PCA, mean_strategies, maxoutdim = 2))
    #println(pca)
    plot!(p, mean_strategies[:, 1], mean_strategies[:, 2], mean_strategies[:, 3])
end
current()

p = plot()
for realization in groupby(adata, :ensemble)
    st = [mean(df.strategy) for df in groupby(realization, :step)]
    mean_strategies = reduce(hcat, st) |> transpose 
    #println(pca)
    plot!(mean_strategies[:, 1])
end
current()



# using ImageCore, CairoMakie
# f = colorsigned() ∘ scalesigned(0., 0, 100.0)
# abm_plot(model, ac = a -> f(a.score))[1]

## not parallel

model = create_model(multiplicative = false, seed = rand(UInt8))

adata, _ = run!(model, player_step!, sampling!, 1000, adata = [:score, :strategy])


result =  groupby(adata, :step)
st = [mean(df.strategy) for df in result]
mean_strategies = reduce(hcat, st) |> transpose 

plot(mean_strategies)

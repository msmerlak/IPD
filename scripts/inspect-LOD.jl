using DrWatson
using DataFramesMeta, StatsBase, StatsPlots
plotly()

include(scriptsdir("run-model.jl"))

@df combine(groupby(run, :step), :fitness => mean) plot(:step, :fitness_mean, xlabel = "Generation", ylabel = "Population mean fitness")
@df combine(groupby(run, :step), :vulnerability => mean) plot!(:step, :vulnerability_mean)

highlander = rand(@subset(run, :step .== T).LOD)
mutations = sign.(reduce(hcat, diff(run[highlander, :].strategy))')

moving_average(vs,n = 10) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]
ma = mapslices(moving_average, mutations; dims = 1)

ma[abs.(ma) .< .8] .= NaN

plot!(abs.(ma))
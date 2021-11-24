

using Distributed
using DataFrames#, Random, CSV, Dates

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


p = Dict(:RSTP => ID_PAYOFFS, :n => n, :t => t, :σ => σ, :multiplicative => multiplicative, :ϵ => 1e-10, :init_share => false, :m => 0.)

# model = create_model(p, rng = MersenneTwister())
models = [create_model(p; rng = MersenneTwister(rand(UInt)), initial_strategy = rand(4)) for _ in 1:PROCESSES] 

adata, mdata = ensemblerun!(models, player_step!, WF_sampling!, 20000, adata = [:strategy, :fitness], parallel = true)

using Query, StatsPlots


plot()
for run in groupby(adata, :ensemble)
    plot!(
        combine(groupby(run, :step), :fitness => mean).fitness_mean
    )
end
current()


run = groupby(adata, :ensemble)[1]
mean_fitness = combine(groupby(run, :step), :fitness => mean).fitness_mean
@gif for t = 1:10:20000
    plts = []
    layout = @layout [grid(2,2)
    b{0.2h}]
    for i in 1:4
        push!(plts, 
        histogram(
            reduce(hcat, run[run.step .== t, :strategy])[i, :],
            legend = false,
            xlabel = "p$i",
            ylabel = "probability",
            xlims = (0, 1),
            ylims = (0, .1),
            bins = 50,
            normalize = :probability,
            color = :black
        )
        )
    end
    push!(plts, 
    plot(mean_fitness[1:t], xlims = (0, 20_000), ylims = (0, 3), line_z = mean_fitness[1:t], legend = false, ylabel = "mean fitness")
    )
    plot(plts..., layout = layout)
end



   
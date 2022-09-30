using DrWatson
@quickactivate "IPD"
include("../src/utils.jl")
include("../src/IPD-ABM.jl")
include("../src/constants.jl")

using StatsPlots

σ, n = 5e-3, 20
init = rand(4)
T = 500_000


p = Dict(:RSTP => [3.0, 0.0, 5.0, 1.0], :n => n, :σ => σ, :init_strategy => init)
model = create_model(p; compute_metrics = false)

run, _ = run!(model, player_step!, WF_sampling!, T, adata = [:strategy, :fitness, :cooperation])

mean_fitness = combine(groupby(run, :step), :fitness => mean).fitness_mean
plot(mean_fitness)

# @df combine(groupby(run, :step), :fitness => mean) plot(:fitness_mean)
# @df combine(groupby(run, :step), :cooperation => mean) plot!(:cooperation_mean)
@gif for t = 6500:1:7000
    plts = []
    push!(plts,
        histogram(
            reduce(hcat, run[run.step .== t, :cooperation])', xlims = (0, 1),
            bins = 20,
            ylims = (0, 20),
            title = "t = $t"),
            )
    plot(plts...)
end


####
@gif for t = 1:100:T
    plts = []
        push!(plts, 
        histogram(
            reduce(hcat, run[run.step .== t, :cooperation])',
            legend = false,
            xlims = (0, 1),
            ylims = (0, 1),
            bins = 20,
            normalize = :probability,
            color = :black
        )
        )
    push!(plts, 
    plot(mean_fitness[1:t], xlims = (0, T), ylims = (0, 3), line_z = mean_fitness[1:t], legend = false, ylabel = "mean fitness")
    )
    plot(plts...)
end


   


####
@gif for t = 1:1000:T
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
            bins = 20,
            normalize = :probability,
            color = :black
        )
        )
    end
    push!(plts, 
    plot(mean_fitness[1:t], xlims = (0, T), ylims = (0, 3), line_z = mean_fitness[1:t], legend = false, ylabel = "mean fitness")
    )
    plot(plts..., layout = layout)
end



   
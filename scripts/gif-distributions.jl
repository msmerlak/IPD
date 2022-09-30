

using StatsPlots

# @df combine(groupby(run, :step), :fitness => mean) plot(:fitness_mean)

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



   
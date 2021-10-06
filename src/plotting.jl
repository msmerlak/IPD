using DrWatson
@quickactivate

using Plots, InteractiveDynamics


function plot_probs(adata; title = nothing)


    p = plot(xlabel = "P(C|CC)", ylabel = "P(C|DC)")
    for realization in groupby(adata, :ensemble)
        st = realization.mean_strategy
        mean_strategies = reduce(hcat, st) |> transpose 
        plot!(p, mean_strategies[:, 1], mean_strategies[:, 3], line_z = adata.step, legend = false)
    end
    p12 = current()

    p = plot(xlabel = "P(C|CD)", ylabel = "P(C|DD)")
    for realization in groupby(adata, :ensemble)
        st = realization.mean_strategy
        mean_strategies = reduce(hcat, st) |> transpose 
        plot!(p, mean_strategies[:, 2], mean_strategies[:, 4], line_z = adata.step, legend = false)
    end
    p34 = current()


    plot(p12, p34, title = title)
    
end
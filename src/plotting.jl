using DrWatson
@quickactivate

using Plots, InteractiveDynamics

function plot_probs(adata)
    result =  groupby(adata, :step)
    n = length(result)

    p = plot(xlabel = "P(C|CC)", ylabel = "P(C|CD)")
    for realization in groupby(adata, :ensemble)
        st = [mean(df.strategy) for df in groupby(realization, :step)]
        mean_strategies = reduce(hcat, st) |> transpose 
        plot!(p, mean_strategies[:, 1], mean_strategies[:, 2], line_z = 1:n)
    end
    p12 = current()

    p = plot(xlabel = "P(C|DC)", ylabel = "P(C|DD)")
    for realization in groupby(adata, :ensemble)
        st = [mean(df.strategy) for df in groupby(realization, :step)]
        mean_strategies = reduce(hcat, st) |> transpose 
        plot!(p, mean_strategies[:, 3], mean_strategies[:, 4], line_z = 1:n)
    end
    p34 = current()


    plot(p12, p34)
    
end


function plot_probs(adata)
    #result =  groupby(adata, :step)
    #n = length(result)

    p = plot(xlabel = "P(C|CC)", ylabel = "P(C|CD)")
    for realization in groupby(adata, :ensemble)
        # st = [df.mean_strategies for df in groupby(realization, :step)]
        st = realization.mean_strategy
        mean_strategies = reduce(hcat, st) |> transpose 
        plot!(p, mean_strategies[:, 1], mean_strategies[:, 2], line_z = adata.step, legend = false)
    end
    p12 = current()

    p = plot(xlabel = "P(C|DC)", ylabel = "P(C|DD)")
    for realization in groupby(adata, :ensemble)
        # st = [df.mean_strategies for df in groupby(realization, :step)]
        st = realization.mean_strategy
        mean_strategies = reduce(hcat, st) |> transpose 
        plot!(p, mean_strategies[:, 3], mean_strategies[:, 4], line_z = adata.step, legend = false)
    end
    p34 = current()


    plot(p12, p34)
    
end
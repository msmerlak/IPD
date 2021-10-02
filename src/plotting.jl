function plot_probs(adata)
    result =  groupby(adata, :step)


    p = plot(xlabel = "P(C|CC)", ylabel = "P(C|CD)")
    for realization in groupby(adata, :ensemble)
        st = [mean(df.strategy) for df in groupby(realization, :step)]
        mean_strategies = reduce(hcat, st) |> transpose 
        plot!(p, mean_strategies[:, 1], mean_strategies[:, 2], color = "green", )
    end
    p12 = current()


    p = plot(xlabel = "P(C|DC)", ylabel = "P(C|DD)")
    for realization in groupby(adata, :ensemble)
        st = [mean(df.strategy) for df in groupby(realization, :step)]
        mean_strategies = reduce(hcat, st) |> transpose 
        plot!(p, mean_strategies[:, 3], mean_strategies[:, 4], color = "green")
    end
    p34 = current()


    plot(p12, p34)
    
end

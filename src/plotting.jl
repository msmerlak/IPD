using DrWatson
@quickactivate

using Plots, CairoMakie, InteractiveDynamics
include(srcdir("memory-one-IPD.jl"))


function plot_strategies_3D(adata, filename)


    p = plot(xlabel = "P(C|CC)", ylabel = "P(C|DC)", zlabel = "Time")
    for d in groupby(adata, :ensemble)
        plot!(
            p,
            [s[1] for s in d.mean_strategy],
            [s[3] for s in d.mean_strategy],
            d.step,
            line_z = d.mean_fitness,
            legend = false,
            dpi = 300,
        )
    end
    p13 = current()

    p = plot(xlabel = "P(C|CD)", ylabel = "P(C|DD)", zlabel = "Time")
    for d in groupby(adata, :ensemble)
        plot!(
            p,
            [s[2] for s in d.mean_strategy],
            [s[4] for s in d.mean_strategy],
            d.step,
            line_z = d.mean_fitness,
            legend = false,
            dpi = 300,
        )
    end
    p24 = current()

    plot(p13, p24)
    savefig(filename)
end

function plot_strategies_2D(adata, filename)


    p = plot(xlabel = "P(C|CC)", ylabel = "P(C|DC)")
    for realization in groupby(adata, :ensemble)
        st = realization.mean_strategy
        mean_strategies = reduce(hcat, st) |> transpose
        plot!(
            p,
            mean_strategies[:, 1],
            mean_strategies[:, 3],
            line_z = adata.step,
            legend = false,
            color = :imola,
            dpi = 300,
        )
    end
    p13 = current()

    p = plot(xlabel = "P(C|CD)", ylabel = "P(C|DD)")
    for realization in groupby(adata, :ensemble)
        st = realization.mean_strategy
        mean_strategies = reduce(hcat, st) |> transpose
        plot!(
            p,
            mean_strategies[:, 2],
            mean_strategies[:, 4],
            line_z = adata.step,
            legend = false,
            color = :imola,
            dpi = 300,
        )
    end
    p24 = current()

    plot(p13, p24)
    savefig(filename)
end


function plot_strategies_1D(adata, filename)

    plots = ()
    labels = ("P(C|CC)", "P(C|CD)", "P(C|DC)", "P(C|DD)")
    for i = 1:4
        p = plot(ylabel = labels[i], xlabel = "Time")
        for d in groupby(adata, :ensemble)
            plot!(
                p,
                d.step,
                [s[i] for s in d.mean_strategy],
                line_z = d.mean_fitness,
                legend = false,
                dpi = 300,
            )
        end
        plots = (plots..., current())
    end

    plt = plot(plots..., layout = 4)
    savefig(filename)
end




function plot_property(adata, property, filename = nothing)

    plt = Plots.plot(ylabel = string(property), xlabel = "Time")
    for d in groupby(adata, :ensemble)
        Plots.plot!(
            plt,
            d.step,
            d[!, property],
            line_z = d.mean_fitness,
            legend = false,
            dpi = 300,
        )
    end
    #plt = current()
    filename !== nothing && savefig(filename)
    current()
end

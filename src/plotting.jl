using DrWatson
@quickactivate

using Plots;gr();
using InteractiveDynamics
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

function plot_std_strategies(adata, filename)



    p = plot(ylabel = "std strategy", xlabel = "Time")
    for d in groupby(adata, :ensemble)
        plot!(
            p,
            d.step,
            d.mean_std_strategy,
            line_z = d.mean_fitness,
            legend = false,
            dpi = 300,
        )
    end



    plt = current()
    savefig(filename)
end



function plot_mean_fitness(adata, filename)



    p = plot(ylabel = "mean fitness", xlabel = "time")
    for d in groupby(adata, :ensemble)
        plot!(
            p,
            d.step,
            d.mean_fitness,
            legend = false,
            dpi = 300,
        )
    end



    plt = current()
    savefig(filename)
end

function plot_extorsion(adata, filename)

    p = plot(ylabel = "extorsion", xlabel = "time")
    for d in groupby(adata, :ensemble)
        plot!(
            p,
            d.step,
            extorsion.(d.mean_strategy),
            legend = false,
            dpi = 300,
        )
    end
    plt = current()
    savefig(filename)
end

function plot_generosity(adata, filename)

    p = plot(ylabel = "generosity", xlabel = "time")
    for d in groupby(adata, :ensemble)
        plot!(
            p,
            d.step,
            generosity.(d.mean_strategy),
            legend = false,
            dpi = 300,
        )
    end
    plt = current()
    savefig(filename)
end
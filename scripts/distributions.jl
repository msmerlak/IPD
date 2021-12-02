
using DrWatson

mkdir(plotsdir("mixed", string(today())))
mkdir(datadir("mixed", string(today())))


begin
    @quickactivate "IPD"
    include("../src/utils.jl")
    include("../src/IPD.jl")
    include("../src/constants.jl")
end

include(srcdir("plotting.jl"))

### distributions

σ, n = 1e-1, 50
init = rand(4)
multiplicative = false


p = Dict(:RSTP => RSTP, :n => n, :t => nothing, :σ => σ, :multiplicative => multiplicative, :ϵ => 1e-10, :init_share => false, :m => 0.)

# model = create_model(p, rng = MersenneTwister())
model = create_model(p)

run, _ = run!(model, player_step!, WF_sampling!, 1_000, adata = [:strategy, :fitness])


mean_fitness = combine(groupby(run, :step), :fitness => mean).fitness_mean
@gif for t = 1:10:1000
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
    plot(mean_fitness[1:t], xlims = (0, 1_000), ylims = (0, 3), line_z = mean_fitness[1:t], legend = false, ylabel = "mean fitness")
    )
    plot(plts..., layout = layout)
end



   
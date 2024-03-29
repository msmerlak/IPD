

using Distributed
using DataFrames, Random, CSV, Dates


mkdir(plotsdir("mixed", string(today())))
mkdir(datadir("mixed", string(today())))

PARALLEL = true
PROCESSES = Sys.CPU_THREADS - 1
PARALLEL && addprocs(PROCESSES)
@everywhere begin
    import Pkg
    Pkg.activate(".")
end
@everywhere using DrWatson
@everywhere begin
    @quickactivate "IPD"
    include("../src/utils.jl")
    include("../src/IPD.jl")
    include("../src/constants.jl")
end

include(srcdir("plotting.jl"))

### different initial conditions
for σ ∈ (1e-1),  n ∈ (100, 500), t ∈ (5, 10, 20), multiplicative ∈ (true, false)

    println("σ = $σ, n = $n, t = $t, multiplicative = $multiplicative")

    p = Dict(:RSTP => RSTP, :n => n, :t => t, :σ => σ, :multiplicative => multiplicative, :ϵ => 1e-s10, :init_share => false, :m => 0.)

    models = [create_model(p; rng = MersenneTwister(rand(UInt)), initial_strategy = rand(4)) for _ in 1:PROCESSES] 

    adata, _ = ensemblerun!(models, player_step!, WF_sampling!, 5000, adata = [(:fitness, mean), (:strategy, mean), (:strategy, mean_std)], 
    parallel = true)


    q = subdict(p, [:n, :t, :σ, :multiplicative])

    CSV.write(datadir("mixed", savename(q, "csv")), adata)

    plot_strategies_1D(adata, plotsdir("mixed", string(today()), savename("1D", q, "png")))

    plot_strategies_2D(adata, plotsdir("mixed", string(today()), savename("2D", q, "png")))

    plot_strategies_3D(adata, plotsdir("mixed", string(today()), savename("3D", q, "png")))
end



### same initial conditions, different seeds (Parallel Evolution, PE)

for σ ∈ (5e-3, 5e-2),  n ∈ (100, 500), t ∈ (5, 10), multiplicative ∈ (true, false)

    println("σ = $σ, n = $n, t = $t, multiplicative = $multiplicative")

    init = rand(4)

    p = Dict(:RSTP => RSTP, :n => n, :t => t, :σ => σ, :multiplicative => multiplicative, :ϵ => 1e-10, :init_share => false, :m => 0.)

    models = [create_model(p; rng = MersenneTwister(x), initial_strategy = init) for x in 1:PROCESSES] 

    adata, _ = ensemblerun!(models, player_step!, WF_sampling!, 5000, adata = [(:fitness, mean), (:strategy, mean), (:strategy, mean_std)], 
    parallel = true)

    q = subdict(p, [:n, :t, :σ, :multiplicative])

    CSV.write(datadir("mixed", string(today()), savename("PE", q, "csv")), adata)

    plot_strategies_1D(adata, plotsdir("mixed", string(today()), savename("PE_1D", q, "png")))

    plot_strategies_2D(adata, plotsdir("mixed", string(today()), savename("PE_2D", q, "png")))

    plot_strategies_3D(adata, plotsdir("mixed", string(today()), savename("PE_3D", q, "png")))
end


### extortionary initial condition 

σ, n, t = 5e-3, 500, 10
init = [11/13, 1/2, 7/26, 0.]
multiplicative = false


p = Dict(:RSTP => RSTP, :n => n, :t => t, :σ => σ, :multiplicative => multiplicative, :ϵ => 1e-10, :init_share => false, :m => 0.)

models = [create_model(p; rng = MersenneTwister(x), initial_strategy = init) for x in 1:PROCESSES] 

adata, _ = ensemblerun!(models, player_step!, WF_sampling!, 20000, adata = [(:fitness, mean), (:strategy, mean), (:strategy, mean_std)], 
parallel = true)

q = subdict(p, [:n, :t, :σ, :multiplicative])

CSV.write(datadir("mixed", string(today()), savename("PE_extortion", q, "csv")), adata)

plot_strategies_1D(adata, plotsdir("mixed", string(today()), savename("PE_extortion_1D", q, "png")))

plot_strategies_2D(adata, plotsdir("mixed", string(today()), savename("PE_extortion_2D", q, "png")))

plot_strategies_3D(adata, plotsdir("mixed", string(today()), savename("PE_extortion_3D", q, "png")))

### mean fitness

σ = 5e-3
n = 500
t = 10
multiplicative = false


println("σ = $σ, n = $n, t = $t, multiplicative = $multiplicative")

p = Dict(:RSTP => RSTP, :n => n, :t => t, :σ => σ, :multiplicative => multiplicative, :ϵ => 1e-10, :init_share => false, :m => 0.)

models = [create_model(p; rng = MersenneTwister(rand(UInt)), initial_strategy = rand(4)) for _ in 1:PROCESSES] 

adata, _ = ensemblerun!(models, player_step!, WF_sampling!, 5000, adata = [(:fitness, mean), (:strategy, mean), (:strategy, mean_std), (:share, mean)], 
parallel = true)

q = subdict(p, [:n, :t, :σ, :multiplicative, :m])


plot_mean_fitness(adata, plotsdir("mixed", string(today()), savename("FITNESS", q, "png")))

### extorsion/generosity

σ = 5e-3
n = 500
t = 10
multiplicative = false

println("σ = $σ, n = $n, t = $t, multiplicative = $multiplicative")

p = Dict(:RSTP => RSTP, :n => n, :t => t, :σ => σ, :multiplicative => multiplicative, :ϵ => 1e-10, :init_share => false, :m => 0.)

models = [create_model(p; rng = MersenneTwister(rand(UInt)), initial_strategy = rand(4)) for _ in 1:PROCESSES] 

adata, _ = ensemblerun!(models, player_step!, WF_sampling!, 10000, adata = [(:fitness, mean), (:strategy, mean), (:cooperation, mean), (:generosity, mean), (:extorsion, mean)], 
parallel = true)

q = subdict(p, [:n, :t, :σ, :multiplicative, :m])

for Y ∈ (:mean_cooperation, :mean_extorsion, :mean_generosity, :mean_fitness)

    plot_property(adata, Y, plotsdir("mixed", string(today()), savename(string(Y), q, "png")))
end





### save whole population

σ = 5e-3
n = 500
t = 10
multiplicative = false


println("σ = $σ, n = $n, t = $t, multiplicative = $multiplicative")

p = Dict(:RSTP => RSTP, :n => n, :t => t, :σ => σ, :multiplicative => multiplicative, :ϵ => 1e-10, :init_share => false, :m => 0.)

models = [create_model(p; rng = MersenneTwister(rand(UInt)), initial_strategy = rand(4)) for _ in 1:PROCESSES] 

adata, _ = ensemblerun!(models, player_step!, WF_sampling!, 5000, adata = [(:fitness, mean), (:strategy, mean_std)], 
parallel = true)

q = subdict(p, [:n, :t, :σ, :multiplicative, :m])

Y = :mean_std_strategy

plot_property(adata, Y, plotsdir("mixed", string(today()), savename(string(Y), q, "png")))

### entropy

σ, n, t = 5e-3, 500, 10
init = rand(4)
multiplicative = false


p = Dict(:RSTP => RSTP, :n => n, :t => t, :σ => σ, :multiplicative => multiplicative, :ϵ => 1e-10, :init_share => false, :m => 0.)

model = create_model(p, rng = MersenneTwister())

adata, mdata = run!(model, player_step!, WF_sampling!, 5000, adata = [:strategy])



using Distributed
using DataFrames, CSV, Dates, Random, RandomNumbers

#PARALLEL = true
PROCESSES = 5#Sys.CPU_THREADS - 1
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
    include("../src/IPD.jl")
    include("../src/utils.jl")
    include("../src/constants.jl")
end

include(srcdir("plotting.jl"))

### distributions

begin
    σ = 1e-2
    n = t = 50
    init = rand(4)
    multiplicative = false
    T = 50_000

    p = Dict(:RSTP => ID_PAYOFFS, :n => n, :t => t, :σ => σ, :multiplicative => multiplicative, :ϵ => 1e-10, :init_share => false, :m => 0.)

    #MersenneTwister(rand(UInt))
    # model = create_model(p, rng = MersenneTwister())
    rng = Random.GLOBAL_RNG
    #rng = RandomNumbers.Xorshifts.Xoroshiro128Star()
    models = [create_model(p; rng = rng, initial_strategy = rand(4)) for _ in 1:PROCESSES] 

    adata, mdata = ensemblerun!(models, player_step!, WF_sampling!, T, adata = [(:fitness, mean), (:strategy, mean)], parallel = true)

end
plot_property(adata, :mean_fitness)

using Revise, DrWatson
@quickactivate "IPD"
includet("../src/utils.jl")
includet("../src/IPD-ABM.jl")
includet("../src/constants.jl")

σ, n = 1e-2, 50
init = WSLS

p = Dict(:RSTP => [3.0, 0.0, 5.0, 1.0], :n => n, :σ => σ, :init_strategy => init)

model = create_model(p; LOD = false, reactive_only = false)
#seed!(model, 11)

T = 100_000
run, _ = run!(model, player_step!, WF_sampling!, T, adata = [:strategy, :fitness, :vulnerability], showprogress=true);
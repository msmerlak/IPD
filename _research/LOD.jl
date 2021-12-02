
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

σ, n = 1e-2, 50
init = rand(4)
multiplicative = false


p = Dict(:RSTP => RSTP, :n => n, :t => nothing, :σ => σ, :multiplicative => multiplicative, :ϵ => 1e-10, :init_share => false, :m => 0., :LOD => true)

# model = create_model(p, rng = MersenneTwister())
model = create_model(p)

run, _ = run!(model, player_step!, WF_sampling!, 1000, adata = [:strategy, :fitness, :LOD])

run[run.step .== 1000, :LOD]

lod = rand(run[run.step .== 10000, :LOD])

str_lod = reduce(hcat, map(id -> run[run.id .== id, :fitness][1], lod))

plot(transpose(str_lod))

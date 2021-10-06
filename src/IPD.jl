using DrWatson
@quickactivate "IPD"

include(srcdir("memory-one-IPD.jl"))

using Agents
using LinearAlgebra:det, dot
using Distributions:Normal, Bernoulli
using Random:MersenneTwister
using NaNMath, StatsBase
#using VSL

const ID_PAYOFFS = [3., 0., 5., 1.]

RAND = [.25, .25, .25, .25]
TFT = [1., 0, 1., 0.]
WSLS = [1., 0., 1., 0.]

mutable struct Mem1Player <: AbstractAgent
    id::Int
    pos::NTuple{2, Int}
    strategy::Vector{Float64}
    scores::Vector{Float64}
    fitness::Float64
    share::Bool
end

function match!((X, Y)::Tuple{Mem1Player, Mem1Player}, model)

    Sx = D(X.strategy, Y.strategy, model.RSTP)/D(X.strategy, Y.strategy, ones(4))
    Sy = D(Y.strategy, X.strategy, model.RSTP)/D(Y.strategy, X.strategy, ones(4))


    if !isfinite(Sx) || !isfinite(Sy)

        Sx = dot(v(X.strategy, Y.strategy), model.RSTP)
        Sy = dot(v(Y.strategy, X.strategy), model.RSTP)

    end

    push!(X.scores, Sx)
    push!(Y.scores, Sy)
    return
end


function create_model(p; initial_strategy = RAND, space = nothing, rng)

    properties = deepcopy(p)

    properties[:init_strategy] = initial_strategy
    if properties[:multiplicative] 
        properties[:selection] = a -> exp(StatsBase.mean(a.scores)) 
    else
        properties[:selection] = a -> StatsBase.mean(a.scores)
    end

    properties[:time] = 0
    
    model = AgentBasedModel(
        Mem1Player, 
        space, 
        properties = properties#, 
        #rng = RandomDevice(properties[:seed])
        #rng =BasicRandomNumberGenerator(VSL_BRNG_MT19937, 12345)
        )

    if model.multiplicative
        @. model.RSTP = log(model.RSTP + properties[:ϵ])
    end

    
    if model.space === nothing 
        for id in 1:model.n
            add_agent!(Mem1Player(id, (1,1), properties[:init_strategy], Float64[1. ], 0., properties[:init_share]), model)
        end
    else
        fill_space!(model, properties[:init_strategy], Float64[1. ], 0., properties[:init_share])
    end
    return model
end

function play_matches!(player, model)
        if model.space === nothing
            competitors = [competitor.second for competitor in rand(model.rng, filter(a -> a.second != player, model.agents), model.t)]
        else
            competitors = nearby_agents(player, model)
        end
        for competitor in competitors
            match!((player, competitor), model)
        end
end

function mutate!(player, model)
    player.strategy = player.strategy + rand(model.rng, Normal(0, model.σ), 4)
    player.strategy = map(s -> window(s, model), player.strategy)
    
    if rand(model.rng, Bernoulli(model.m))
        player.share = !player.share
    end
    return
end

function window(x, model)
    if 0. <= x <= 1.
        return x
    elseif x < 0.
        return 1e-6*rand(model.rng)
    elseif x > 1.
        return 1. - 1e-6*rand(model.rng)
    end
end

function player_step!(player, model)
    mutate!(player, model)
    play_matches!(player, model)
end


#include(srcdir("sampling.jl"))
function WF_sampling!(model)

    for a in values(model.agents)
        a.fitness = model.selection(a)
        a.scores = Float64[]
    end



    Agents.sample!(model, model.n, :fitness)
    model.time += 1
end
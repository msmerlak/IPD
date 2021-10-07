using DrWatson
@quickactivate "IPD"

include(srcdir("memory-one-IPD.jl"))

using Agents
using Distributions:Normal, Bernoulli
using Random:MersenneTwister
using NaNMath, StatsBase
#using VSL



mutable struct Mem1Player <: AbstractAgent
    id::Int
    pos::NTuple{2, Int}
    strategy::Vector{Float64}
    scores::Vector{Float64}
    fitness::Float64
    share::Bool
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
        properties = properties, 
        rng = rng
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


function match!((X, Y)::Tuple{Mem1Player, Mem1Player}, model)
    
    Z = Δ(X.strategy, Y.strategy, ones(4))
    Sx = Δ(X.strategy, Y.strategy, model.RSTP)/Z
    Sy = Δ(Y.strategy, X.strategy, model.RSTP)/Z


    if !isfinite(Sx) || !isfinite(Sy)

        Sx = dot(v(X.strategy, Y.strategy), model.RSTP)
        Sy = dot(v(Y.strategy, X.strategy), model.RSTP)

    end

    push!(X.scores, Sx)
    push!(Y.scores, Sy)
    return 
end


function play_matches!(player, model)
        if model.space === nothing
            competitors = [competitor.second for competitor in filter(a -> a.second != player, rand(model.rng, model.agents, model.t))]
        else
            competitors = nearby_agents(player, model)
        end
        
        for competitor in competitors
            match!((player, competitor), model)
        end
end


function mutate!(player, model)

    player.strategy += rand(model.rng, Normal(0, model.σ), 4)
    window!(player.strategy, model)

    if model.m > 0 && rand(model.rng, Bernoulli(model.m))
        player.share = !player.share
    end
    return
end

function window!(x, model)
        @. x = min(max(x, 1e-6*rand(model.rng)), 1. - 1e-6*rand(model.rng))
    return
end

function player_step!(player, model)
    mutate!(player, model)
    play_matches!(player, model)
end


#include(srcdir("sampling.jl"))
function WF_sampling!(model)

    # compute fitness
    for a in allagents(model)
        a.fitness = model.selection(a)
        a.scores = Float64[]
    end

    # share fitness among sharers
    sharers = (a for a in allagents(model) if a.share)
    non_sharers = (a for a in allagents(model) if !a.share)
    pooled_fitness = mean([a.fitness for a in sharers])
    for a in sharers
        a.fitness = pooled_fitness
    end

   

    Agents.sample!(model, model.n, :fitness)
    model.time += 1

end
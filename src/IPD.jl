using DrWatson
@quickactivate "IPD"

include(srcdir("memory-one-IPD.jl"))

using Agents
using LinearAlgebra:det, dot
using Distributions:Normal
using Random:MersenneTwister
using KrylovKit:eigsolve
using StatsBase:mean

const ID_PAYOFFS = [3., 0., 5., 1.]

RAND = [.25, .25, .25, .25]
TFT = [1., 0, 1., 0.]


mutable struct Mem1Player <: AbstractAgent
    id::Int
    pos::NTuple{2, Int}
    strategy::Vector{Float64}
    score::Vector{Float64}
end

function match!((X, Y)::Tuple{Mem1Player, Mem1Player}, model)

    if X.strategy != Y.strategy
        Z = D(X.strategy, Y.strategy, ones(4))
        push!(X.score, D(X.strategy, Y.strategy, model.RSTP)/Z)
        push!(Y.score, D(Y.strategy, X.strategy, model.RSTP)/Z)
    else # other way of computing the same score that also works when identical strategies meet
        push!(X.score, dot(v(X.strategy, Y.strategy), model.RSTP))
        push!(Y.score, dot(v(Y.strategy, X.strategy), model.RSTP))
    end

    return
end


function create_model(p; space = nothing)

    properties = deepcopy(p)

    if properties[:multiplicative] 
        properties[:selection] = a -> exp(mean(a.score))
    else
        properties[:selection] = a -> mean(a.score)
    end

    
    model = AgentBasedModel(
        Mem1Player, 
        space, 
        properties = properties, 
        rng = MersenneTwister(properties[:seed]))

    if model.multiplicative
        @. model.RSTP = log(model.RSTP + properties[:ϵ])
    end

    
    if model.space === nothing 
        for id in 1:model.n
            add_agent!(Mem1Player(id, (1,1), properties[:init], Float64[]), model)
        end
    else
        fill_space!(model, properties[:init], Float64[])
    end
    return model
end

function play_matches!(player, model)
        if model.space === nothing
            competitors = [competitor.second for competitor in rand(model.rng, model.agents, model.t)]
        else
            competitors = nearby_agents(player, model)
        end
        for competitor in competitors
            match!((player, competitor), model)
        end
end

function mutate!(player, model)
    player.strategy = player.strategy + rand(model.rng, Normal(0, model.σ), 4)
    player.strategy = window.(player.strategy)
    return
end

function window(x)
    return min(max(x, 0.), 1.)
end

function player_step!(player, model)
    play_matches!(player, model)
    mutate!(player, model)
end

function WF_sampling!(model)
    # Wright-Fisher sampling
    Agents.sample!(model, model.n, model.selection)

    # Reset scores
    for a in values(model.agents)
        a.score = Float64[]
    end
end
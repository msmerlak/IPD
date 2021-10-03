using DrWatson
@quickactivate "IPD"

include("memory-one-IPD.jl")

using Agents
using LinearAlgebra:det, dot
using Distributions:Normal
using Random:MersenneTwister
using KrylovKit:eigsolve
using StatsBase:mean

RAND = [.25, .25, .25, .25]
TFT = [1., 0, 1., 0.]

mutable struct Mem1Player <: AbstractAgent
    id::Int
    pos::NTuple{2, Int}
    strategy::Vector{Float64}
    score::Vector{Float64}
end

function match!(players::Tuple{Mem1Player, Mem1Player}, model)
    X, Y = players
    p, q = X.strategy, Y.strategy
    R, S, T, P = model.RSTP
    S₁= [R, S, T, P]
    S₂= [R, T, S, P]

    if p != q
        push!(X.score, D(p, q, S₁)/D(p, q, ones(length(p))))
        push!(Y.score, D(p, q, S₂)/D(p, q, ones(length(p))))
    else # other way of computing the same score that also works when identical strategies meet
        push!(X.score, dot(v(p, q), S₁))
        push!(Y.score, dot(v(p, q), S₂))
    end

    return
end

function create_model(; 
    multiplicative = false, 
    ϵ = 1e-10, 
    RSTP = [3., 0., 5., 1.], 
    popsize = 500, 
    tournament_size = 10, 
    mutational_effect = 1e-2, 
    space = nothing, 
    initial_strategy = RAND, 
    reproduction_rate = 1e-1,
    seed = 1)

    if multiplicative 
        selection_function = a -> exp(mean(a.score))
    else
        selection_function = a -> mean(a.score)
    end

    properties = Dict(
        :RSTP => RSTP, 
        :n => popsize, 
        :t => tournament_size, 
        :σ => mutational_effect, 
        :multiplicative => multiplicative, 
        :init => initial_strategy, 
        :r => reproduction_rate,
        :selection => selection_function)
    
    model = AgentBasedModel(
        Mem1Player, 
        space, 
        properties = properties, 
        rng = MersenneTwister(seed))

    if multiplicative
        @. RSTP = log(RSTP + ϵ)
    end

    
    
    if space === nothing 
        for id in 1:popsize
            add_agent!(Mem1Player(id, (1,1), initial_strategy, Float64[]), model)
        end
    else
        fill_space!(model, initial_strategy, Float64[])
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
    return min(max(x, 1e-6), 1-1e-6)
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
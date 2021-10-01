using DrWatson
@quickactivate "IPD"

using Agents
using LinearAlgebra:det
using Distributions:Normal
using Random:MersenneTwister
using KrylovKit:eigsolve

const RAND = [.25, .25, .25, .25]
const TFT = [1., 0, 1., 0.]

mutable struct Mem1Player <: AbstractAgent
    id::Int
    pos::NTuple{2, Int}
    strategy::Vector{Float64}
    score::Float64
end

## Press-Dyson determinant
function D(p, q, f)
    return @fastmath det([[-1 + p[1]*q[1], - 1 + p[1], - 1 + q[1], f[1]] [p[2]*q[3], -1 + p[2], q[3], f[2]] [p[3]*q[2], p[3], -1+q[2], f[3]] [p[4]*q[4], p[4], q[4], f[4]]])
end

function M(p, q)
    qq = [q[1], q[3], q[2], q[4]]
    return @. [p*qq p*(1-qq) (1-p)*qq (1-p)*(1-qq)]
end

function v(p, q)
    val, vec, conv  = eigsolve(transpose(M(p, q)), 1, :LR)
    return vec[1]./sum(vec[1])
end

function match!(players::Tuple{Mem1Player, Mem1Player}, model)
    X, Y = players
    p, q = X.strategy, Y.strategy
    R, S, T, P = model.RSTP
    S₁= [R, S, T, P]
    S₂= [R, T, S, P]
    #println(D(p, q, S₁)/D(p, q, ones(length(p))))
    #M = @. [p*q p*(1-q) (1-p)*q (1-p)*(1-q)] # Markov matrix for mem-1 IPD
    #val, vec, conv  = eigsolve(transpose(M), 1, :LR)
    #v = vec[1]
    X.score += D(p, q, S₁)/D(p, q, ones(length(p)))
    Y.score += D(p, q, S₂)/D(p, q, ones(length(p)))

    return
end

function create_model(; multiplicative = false, ϵ = 0., RSTP = [3., 0., 5., 1.], popsize = 100, tournament_size = 10, mutational_effect = 1e-2, space = nothing, seed = 1, initial_strategy = RAND)


    properties = Dict(:RSTP => RSTP, :n => popsize, :t => tournament_size, :σ => mutational_effect, :multiplicative => multiplicative, :init => initial_strategy)
    
    model = AgentBasedModel(Mem1Player, space, properties = properties, rng = MersenneTwister(seed))

    if multiplicative
        @. RSTP = log(RSTP + ϵ)
    end

    initial_score = 0.
    
    if space === nothing 
        for id in 1:popsize
            add_agent!(Mem1Player(id, (1,1), initial_strategy, initial_score), model)
        end
    else
        fill_space!(model, initial_strategy, initial_score)
    end
    return model
end

function play_matches!(player, model)
        if model.space === nothing
            for competitor in rand(model.rng, model.agents, model.t)
                match!((player, competitor.second), model)
            end
        else
            for competitor in nearby_agents(player, model)
                match!((player, competitor), model)
            end
        end
end

function mutate!(player, model)
    player.strategy = player.strategy + rand(model.rng, Normal(0, model.σ), 4)
    player.strategy = window.(player.strategy)
    return
end

function window(x)
    return min(max(x, 0), 1)
end


function player_step!(player, model)
    player.score = 0.
    play_matches!(player, model)
    if model.multiplicative
        player.score = exp(player.score)
    end
    mutate!(player, model)
    if model.space !== nothing
        walk!(player, rand, model)
    end
end

function WF_sampling!(model)
    Agents.sample!(model, model.n, :score)
end



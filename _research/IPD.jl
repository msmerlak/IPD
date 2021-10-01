using DrWatson
@quickactivate "IPD"

using Agents
using LinearAlgebra:det, dot
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

# p = P(C|CC, CD, DC, DD)

function D(p, q, f)
    return @fastmath det([[-1 + p[1]*q[1], - 1 + p[1], - 1 + q[1], f[1]] [p[2]*q[3], -1 + p[2], q[3], f[2]] [p[3]*q[2], p[3], -1+q[2], f[3]] [p[4]*q[4], p[4], q[4], f[4]]])
end

function M(p, q)
    qq = [q[1], q[3], q[2], q[4]]
    return @. [p*qq p*(1-qq) (1-p)*qq (1-p)*(1-qq)]
end

function v(p, q)
    val, vec, conv  = eigsolve(transpose(M(p, q)), 1, :LR)
    return real(vec[1]./sum(vec[1]))
end

function match!(players::Tuple{Mem1Player, Mem1Player}, model)
    X, Y = players
    p, q = X.strategy, Y.strategy
    R, S, T, P = model.RSTP
    S₁= [R, S, T, P]
    S₂= [R, T, S, P]

    if p != q
        X.score += D(p, q, S₁)/D(p, q, ones(length(p)))
        Y.score += D(p, q, S₂)/D(p, q, ones(length(p)))
    else # other way of computing the same score that also works with degenerate strategies
        X.score += dot(v(p, q), S₁)
        Y.score += dot(v(p, q), S₂)
    end

    # if D(p, q, S₁) > D(p, q, S₂) && rand(model.rng) < model.r
    #     add_agent!(model, X.pos, X.strategy, X.score)
    # end
    return
end

function create_model(; multiplicative = false, ϵ = 0., RSTP = [3., 0., 5., 1.], popsize = 500, tournament_size = 10, mutational_effect = 1e-2, space = nothing, seed = 1, initial_strategy = RAND, reproduction_rate = 1e-1)


    properties = Dict(:RSTP => RSTP, :n => popsize, :t => tournament_size, :σ => mutational_effect, :multiplicative => multiplicative, :init => initial_strategy, :r => reproduction_rate)
    
    model = AgentBasedModel(Mem1Player, space, properties = properties, rng = MersenneTwister(seed))

    if multiplicative
        @. RSTP = log(RSTP + ϵ)
    end

    initial_score = .1
    
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

    player.score = 1e-6
    play_matches!(player, model)
    if model.multiplicative
        player.score = exp(player.score)
    end

    mutate!(player, model)
    
    # if model.space !== nothing
    #     walk!(player, rand, model)
    # end
end

function sampling!(model)
    #Agents.sample!(model, model.n)
    Agents.sample!(model, model.n, a -> exp(a.score))
end
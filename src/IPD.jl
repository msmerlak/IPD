using DrWatson
@quickactivate "IPD"

include(srcdir("memory-one-IPD.jl"))

using Agents
using Distributions: Normal, Bernoulli
using Random: MersenneTwister
using StatsBase: mean
import Entropies

mutable struct Mem1Player <: AbstractAgent
    id::Int
    pos::NTuple{2,Int}
    strategy::Vector{Float64}
    scores::Vector{Float64}
    fitness::Float64
    share::Bool
    cooperation::Float64
    generosity::Float64
    extorsion::Float64
end

function create_model(
    p;
    initial_strategy = RAND,
    space = nothing,
    compute_metrics = false,
    rng = Random.GLOBAL_RNG,
)

    properties = deepcopy(p)

    properties[:init_strategy] = initial_strategy
    properties[:compute_metrics] = compute_metrics

    if properties[:multiplicative]
        properties[:selection] = a -> exp(mean(a.scores))
    else
        properties[:selection] = a -> mean(a.scores)
    end

    model = AgentBasedModel(Mem1Player, space, properties = properties, rng = rng)

    if model.multiplicative
        @. model.RSTP = log(model.RSTP + properties[:ϵ])
    end


    if model.space === nothing
        for id = 1:model.n
            add_agent!(
                Mem1Player(
                    id,
                    (1, 1),
                    properties[:init_strategy],
                    Float64[1.0],
                    0.0,
                    properties[:init_share],
                    0.0,
                    0.0,
                    0.0,
                ),
                model,
            )
        end
    else
        fill_space!(
            model,
            properties[:init_strategy],
            Float64[1.0],
            0.0,
            properties[:init_share],
            0.0,
            0.0,
            0.0,
        )
    end
    return model
end


function match!((X, Y)::Tuple{Mem1Player,Mem1Player}, model)

    ## Press-Dyson formula
    #if X.strategy != Y.strategy
        push!(X.scores, π(X.strategy, Y.strategy, model.RSTP))
        push!(Y.scores, π(Y.strategy, X.strategy, model.RSTP))
    #end
end


function play_matches!(player, model)

    ## pick competitors
    if model.space === nothing
        competitors = allagents(model)
    else
        competitors = nearby_agents(player, model)
    end

    ## play matches
    for competitor in competitors
        match!((player, competitor), model)
    end

    ## compute metrics
    #player.cooperation = cooperation(player, competitors)
    #player.generosity = generosity(player, competitors)
    #player.extorsion = extorsion(player, competitors)
end


function mutate!(player, model)

    player.strategy .+= rand(model.rng, Normal(0., model.σ), 4)
    
    window!(player.strategy, model)
    #player.total_mutations += sum(new_strategy - player.strategy)


    # if model.m > 0 && rand(model.rng, Bernoulli(model.m))
    #     player.share = !player.share
    # end

end

function window!(x, model)
    @. x = min(max(x, 0.), 1.0)
    #@. x = min(max(x, 1e-6 * rand(model.rng)), 1.0 - 1e-6 * rand(model.rng))
    return
end

function player_step!(player, model)
    mutate!(player, model)
    play_matches!(player, model)
end



#include(srcdir("sampling.jl"))
function WF_sampling!(model)

    
    # compute fitness and reset scores etc
    for a in allagents(model)
        a.fitness = model.selection(a)
        a.scores = Float64[]
        sizehint!(a.scores, model.n)
        #a.outgoing_mutations = 0
    end

    # share fitness among sharers
    # sharers = (a for a in allagents(model) if a.share)
    # non_sharers = (a for a in allagents(model) if !a.share)
    # pooled_fitness = mean([a.fitness for a in sharers])
    # for a in sharers
    #     a.fitness = pooled_fitness
    # end


    Agents.sample!(model, model.n, :fitness)

end

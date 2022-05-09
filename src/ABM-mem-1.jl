using DrWatson

include(srcdir("memory-one-IPD.jl"))
include(srcdir("sample-with-LOD.jl"))
include(srcdir("constants.jl"))

using Agents
using Distributions: Normal, Bernoulli
using Random: GLOBAL_RNG
using StatsBase: mean
using StaticArrays


mutable struct Mem1Player <: AbstractAgent
    id::Int
    pos::NTuple{2,Int}
    strategy::MVector{4, Float64}
    scores::Vector{Float64}
    fitness::Float64
    LOD::Vector{Int64}
    cooperation::Float64
    generosity::Float64
    extorsion::Float64
end

function create_model(
    p;
    space = nothing,
    compute_metrics = false,
    LOD = false,
    RSTP = RSTP,
    rng = GLOBAL_RNG
)

    properties = deepcopy(p)

    properties[:compute_metrics] = compute_metrics
    properties[:LOD] = LOD
    properties[:RSTP] = RSTP
    properties[:selection] = a -> mean(a.scores)

    model = AgentBasedModel(Mem1Player, space, properties = properties, rng = rng)
    model.n = Int(model.n)


    if model.space === nothing
        for id = 1:model.n
            add_agent!(
                Mem1Player(
                    id,
                    (1, 1),
                    model.init_strategy,
                    Float64[1.0],
                    0.0,
                    Int64[],
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
            model.init_strategy,
            Float64[1.0],
            0.0,
            Int64[],
            0.0,
            0.0,
            0.0,
        )
    end
    return model
end


function match!((X, Y)::Tuple{Mem1Player,Mem1Player}, model)

    ## Press-Dyson formula for stationary payoffs
    if X.id != Y.id
        push!(X.scores, π(X.strategy, Y.strategy, model.RSTP))
        push!(Y.scores, π(Y.strategy, X.strategy, model.RSTP))
    end
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
    if model.compute_metrics
        player.cooperation = cooperation(player, competitors)
        #player.generosity = generosity(player, competitors)
        #player.extorsion = extorsion(player, competitors)
    end
end


function mutate!(player, model)
    jiggle!(player.strategy, model.σ, model.rng)
    window!(player.strategy)
end

function jiggle!(x::MVector{4}, σ, rng)
    x .+= σ * (2*rand(rng, MVector{4}) .- 1)
end

function window!(x)
    @. x = min(max(x, 1e-9), 1.0 - 1e-9)
end

function player_step!(player, model)
    mutate!(player, model)
    play_matches!(player, model)
end

function WF_sampling!(model)

    # compute fitness and reset scores etc
    for a in allagents(model)
        a.fitness = model.selection(a)
        a.scores = Float64[]
        sizehint!(a.scores, model.n^2)
    end

    # Wright-Fisher sampling
    if !model.LOD
        Agents.sample!(model, model.n, :fitness)
    else
        sample_with_LOD!(model, model.n, :fitness)
    end

end

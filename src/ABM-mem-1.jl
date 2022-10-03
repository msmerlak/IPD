using DrWatson

include(srcdir("Mem1-IPD.jl"))
include(srcdir("sample-with-LOD.jl"))
include(srcdir("constants.jl"))

using Agents
using Random: GLOBAL_RNG
using StatsBase: mean
using StaticArrays

mutable struct Mem1Player <: AbstractAgent
    id::Int
    pos::NTuple{2,Int}
    strategy::MVector{4,Float64}
    scores::Vector{Float64}
    fitness::Float64
    LOD::Vector{Int64}
    vulnerability::Float64
end

function create_model(
    p;
    space=nothing,
    LOD=false,
    reactive_only=false
)

    properties = deepcopy(p)

    properties[:reactive_only] = reactive_only
    properties[:LOD] = LOD
    properties[:fitness] = a -> mean(a.scores)

    model = AgentBasedModel(
        Mem1Player, space;
        properties=properties)
    model.n = Int(model.n)

    if model.space === nothing
        for id = 1:model.n
            add_agent!(
                Mem1Player(
                    id,
                    (1, 1),
                    model.init_strategy,
                    Float64[1.0],
                    NaN,
                    Int64[],
                    NaN
                ),
                model,
            )
        end
    else
        fill_space!(
            model,
            model.init_strategy,
            Float64[1.0],
            NaN,
            Int64[],
            NaN
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
end


function mutate!(player, model)
    player.strategy .+= σ * (2 * rand(model.rng, MVector{4}) .- 1)
    chop!(player.strategy)
end


function chop!(x, ϵ=1e-9)
    @. x = min(max(x, ϵ), 1.0 - ϵ)
end

function player_step!(player, model)
    mutate!(player, model)
    player.vulnerability = 1 - robustness(player, model)
    if model.reactive_only
        player.strategy[[1, 3]] .= player.strategy[[2, 4]]
    end
    play_matches!(player, model)
end

function WF_sampling!(model)

    # compute fitness and reset scores etc
    for a in allagents(model)
        a.fitness = model.fitness(a)
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

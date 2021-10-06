using Agents
using Test
using Random

mutable struct Agent1 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent2 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent3 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent4 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent5 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent6 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent7 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent8 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent9 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent10 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent11 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent12 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent13 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent14 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct Agent15 <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    money::Int
end

mutable struct BadAgent <: AbstractAgent
    useless::Int
    pos::Int
end
mutable struct BadAgentId <: AbstractAgent
    id::Float64
end
struct ImmutableAgent <: AbstractAgent
    id::Int
end
mutable struct DiscreteVelocity <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Int}
    diameter::Float64
end
mutable struct ParametricAgent{T<:Integer} <: AbstractAgent
    id::T
    pos::NTuple{2,T}
    weight::T
    info::String
end

@testset "sample!" begin
    Random.seed!(6459)
    model = ABM(Agent2)
    for i in 1:20
        add_agent!(model, rand(model.rng) / rand(model.rng))
    end
    allweights = [i.weight for i in values(model.agents)]
    mean_weights = sum(allweights) / length(allweights)
    sample!(model, 12, :weight)
    @test nagents(model) == 12
    allweights = [i.weight for i in values(model.agents)]
    mean_weights_new = sum(allweights) / length(allweights)
    @test mean_weights_new > mean_weights
    sample!(model, 40, :weight)
    @test nagents(model) == 40
    allweights = [i.weight for i in values(model.agents)]
    mean_weights_new = sum(allweights) / length(allweights)
    @test mean_weights_new > mean_weights

    Random.seed!(6459)
    model2 = ABM(Agent3, GridSpace((10, 10)))
    for i in 1:20
        add_agent_single!(Agent3(i, (1, 1), rand(model2.rng) / rand(model2.rng)), model2)
    end
    @test sample!(model2, 10) === nothing
    @test sample!(model2, 10, :weight) === nothing
    allweights = [i.weight for i in values(model2.agents)]
    mean_weights = sum(allweights) / length(allweights)
    sample!(model2, 12, :weight)
    @test nagents(model2) == 12
    allweights = [i.weight for i in values(model2.agents)]
    mean_weights_new = sum(allweights) / length(allweights)
    @test mean_weights_new > mean_weights

    sample!(model2, 40, :weight)
    @test nagents(model2) == 40

    Random.seed!(6459)
    #Guarantee all starting weights are unique
    model3 = ABM(Agent2)
    while true
        for i in 1:20
            add_agent!(model3, rand(model3.rng) / rand(model3.rng))
        end
        allweights = [i.weight for i in values(model3.agents)]
        allunique(allweights) && break
    end
    # Cannot draw 50 samples out of a pool of 20 without replacement
    @test_throws ErrorException sample!(model3, 50, :weight; replace = false)
    sample!(model3, 15, :weight; replace = false)
    allweights = [i.weight for i in values(model3.agents)]
    @test allunique(allweights)
    model3 = ABM(Agent2)
    while true
        for i in 1:20
            add_agent!(model3, rand(model3.rng) / rand(model3.rng))
        end
        allweights = [i.weight for i in values(model3.agents)]
        allunique(allweights) && break
    end
    sample!(model3, 100, :weight; replace = true)
    allweights = [i.weight for i in values(model3.agents)]
    @test !allunique(allweights)
end
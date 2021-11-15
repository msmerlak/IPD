using DrWatson
import Statistics:std

#import Clustering:dbscan
using Agents

function subdict(dict, keys)
    subdict = Dict()
    for (k, v) in dict
        k in keys && push!(subdict, k => v)
    end
    return subdict
end

function std(X::T) where T<: Base.Generator
    std([x for x in X])
end

function mean_std(X::T) where T<: Base.Generator
    mean(std([x for x in X]))
end

# function entropy(model::AgentBasedModel)
#     genentropy(Dataset([player.strategy for player in allagents(model)]), NaiveKernel(.1)
#     )
# end

# function density(model::AgentBasedModel)
#     probabilities(Dataset([player.strategy for player in allagents(model)]), 2
#     )
# end


# function clust(model)
#     clstr = dbscan(reduce(hcat, [player.strategy for player in allagents(model)]), .05, min_cluster_size = 10)
#     return clstr#[c.size for c in clstr]
# end
using DrWatson
@quickactivate "IPD"

import Statistics:std, var
#import Clustering:dbscan
using Agents


#mean(x::Vector) = mean(filter(!isnan, x))


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


function mean_cooperation(adata::DataFrame)
    mean(adata.mean_cooperation)
end


function var_cooperation(adata::DataFrame)
    mean(adata.cooperation)
end
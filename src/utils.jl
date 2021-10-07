using DrWatson
import Statistics:std

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
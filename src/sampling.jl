### Adapted from Agents.sample! which is buggy

using Distributions:wsample
using Agents:add_newids!


get_data(a, s::Symbol, obtainer::Function = identity) = obtainer(getproperty(a, s))
get_data(a, f::Function, obtainer::Function = identity) = obtainer(f(a))


function mysample!(
    model::ABM,
    n::Int,
    weight = nothing;
    replace = true,
)
    nagents(model) > 0 || return
    org_ids = collect(keys(model.agents))
    if weight !== nothing
        weights = [get_data(a, weight, identity) for a in values(model.agents)]
        newids = wsample(model.rng, org_ids, weights, n, replace = replace)
    else
        newids = wsample(model.rng, org_ids, n, replace = replace)
    end
    add_newids!(model, org_ids, newids)
end
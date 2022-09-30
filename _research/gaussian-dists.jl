using Distributions, Expectations
using Cubature


function Π(p::Vector{Float64}, q::Vector{Float64}, σ::Float64)
    int, err = hcubature(x -> π(x[1:4], x[5:8])*exp(-sum((x[1:4] - p).^2 + (x[5:8] - q).^2)/(2σ^2)), zeros(8), ones(8); reltol = 1e-3, abstol = 1e-3)
    return int*(2*Base.π*σ^2)^(-4), err
end


p, q = ones(4), ones(4)
Π(q, q, .2)
π(p, q)
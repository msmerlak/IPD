using DrWatson

using KrylovKit:eigsolve
using LinearAlgebra:det, dot
include(srcdir("constants.jl"))
## Press-Dyson determinant

# p = P(C|CC, CD, DC, DD)

function D(p::Vector{Float64}, q::Vector{Float64}, f::Vector{Float64})
    return @fastmath det([[-1 + p[1]*q[1], - 1 + p[1], - 1 + q[1], f[1]] [p[2]*q[3], -1 + p[2], q[3], f[2]] [p[3]*q[2], p[3], -1+q[2], f[3]] [p[4]*q[4], p[4], q[4], f[4]]])
end

function M(p, q)
    qq = [q[1], q[3], q[2], q[4]]
    return @. [p*qq p*(1-qq) (1-p)*qq (1-p)*(1-qq)]
end

function v(p::Vector{Float64}, q::Vector{Float64})
    val, vec, conv  = eigsolve(transpose(M(p, q)), 1, :LR)
    return real(vec[1]./sum(vec[1]))
end

## faster, using exact result

function Δ(p::Vector{Float64}, q::Vector{Float64}, f::Vector{Float64})
    return (p[1]*q[1] - 1)*((f[3]*q[3] - (q[2] - 1)*f[2])*p[4] + (p[2] - 1)*((q[2] - 1)*f[4] - f[3]*q[4]) - (f[4]*q[3] - f[2]*q[4])*p[3]) + ((1 - p[2])*((q[1] - 1)*f[4] - f[1]*q[4]) + ((q[1] - 1)*f[2] - f[1]*q[3])*p[4] + (p[1] - 1)*(f[4]*q[3] - f[2]*q[4]))*p[3]*q[2] - ((1 - p[2])*((q[1] - 1)*f[3] - (q[2] - 1)*f[1]) + ((q[1] - 1)*f[2] - f[1]*q[3])*p[3] + (p[1] - 1)*(f[3]*q[3] - (q[2] - 1)*f[2]))*p[4]*q[4] - ((p[1] - 1)*((q[2] - 1)*f[4] - f[3]*q[4]) + ((q[1] - 1)*f[3] - (q[2] - 1)*f[1])*p[4] - ((q[1] - 1)*f[4] - f[1]*q[4])*p[3])*p[2]*q[3]
end

function π(p::Vector{Float64}, q::Vector{Float64})
    return Δ(p, q, ID_PAYOFFS)/Δ(p, q, ones(4))
end

R, S, T, P = ID_PAYOFFS
E(p, q) = (π(p,q) - P)/(π(q,p) - P)
G(p, q) = (R-π(p,q))/(R - π(q,p)) 

function absolute_extorsion(p::Vector{Float64})
    χ = map(q-> E(p,q) - 1, [rand(4) for _ in 1:100])
    return max(minimum(χ), 0)
end

function absolute_generosity(p::Vector{Float64})
    χ = map(q-> G(p,q) - 1, [rand(4) for _ in 1:100])
    return max(minimum(χ), 0)
end

function extorsion(player, competitors::Base.Generator)
    χ = map(q-> E(player.strategy,q) - 1, [c.strategy for c in competitors])
    return max(minimum(χ), 0)
end

function generosity(player, competitors::Base.Generator)
    χ = map(q-> G(player.strategy,q) - 1, [c.strategy for c in competitors])
    return max(minimum(χ), 0)
end

function cooperation(player, competitors::Base.Generator)
    χ = map(q-> Δ(player.strategy, q, [1., 1., 0., 0.])/Δ(player.strategy ,q, ones(4)), [c.strategy for c in competitors])
    return mean(χ)
end

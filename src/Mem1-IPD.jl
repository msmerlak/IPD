using DrWatson

#include(srcdir("IPD.jl"))

#using KrylovKit: eigsolve
using LinearAlgebra: det, dot
using DataFrames

include(srcdir("constants.jl"))

## Press-Dyson determinant
# p = P(C|CC, CD, DC, DD)

function D(p::Vector{Float64}, q::Vector{Float64}, f::Vector{Float64})
    return @fastmath det(
        [[-1 + p[1] * q[1], -1 + p[1], -1 + q[1], f[1]] [
            p[2] * q[3],
            -1 + p[2],
            q[3],
            f[2],
        ] [p[3] * q[2], p[3], -1 + q[2], f[3]] [p[4] * q[4], p[4], q[4], f[4]]],
    )
end


function M(p, q)
    qq = [q[1], q[3], q[2], q[4]]
    return @. [p * qq p * (1 - qq) (1 - p) * qq (1 - p) * (1 - qq)]
end

function v(p::Vector{Float64}, q::Vector{Float64})
    val, vec, conv = eigsolve(transpose(M(p, q)), 1, :LR)
    return real(vec[1] ./ sum(vec[1]))
end

## faster, using exact result

# function Δ(p::Vector, q::Vector, f::Vector)
#     return (p[1] * q[1] - 1) * (
#         (f[3] * q[3] - (q[2] - 1) * f[2]) * p[4] +
#         (p[2] - 1) * ((q[2] - 1) * f[4] - f[3] * q[4]) - (f[4] * q[3] - f[2] * q[4]) * p[3]
#     ) +
#            (
#                (1 - p[2]) * ((q[1] - 1) * f[4] - f[1] * q[4]) +
#                ((q[1] - 1) * f[2] - f[1] * q[3]) * p[4] +
#                (p[1] - 1) * (f[4] * q[3] - f[2] * q[4])
#            ) *
#            p[3] *
#            q[2] -
#            (
#                (1 - p[2]) * ((q[1] - 1) * f[3] - (q[2] - 1) * f[1]) +
#                ((q[1] - 1) * f[2] - f[1] * q[3]) * p[3] +
#                (p[1] - 1) * (f[3] * q[3] - (q[2] - 1) * f[2])
#            ) *
#            p[4] *
#            q[4] -
#            (
#                (p[1] - 1) * ((q[2] - 1) * f[4] - f[3] * q[4]) +
#                ((q[1] - 1) * f[3] - (q[2] - 1) * f[1]) * p[4] -
#                ((q[1] - 1) * f[4] - f[1] * q[4]) * p[3]
#            ) *
#            p[2] *
#            q[3]
# end
# function  π(p, q)
#     return Δ(p, q, RSTP) / Δ(p, q, ones(4))

# end

function π(p, q, f = RSTP)
    return ((p[1] * q[1] - 1) * (
        (f[3] * q[3] - (q[2] - 1) * f[2]) * p[4] +
        (p[2] - 1) * ((q[2] - 1) * f[4] - f[3] * q[4]) - (f[4] * q[3] - f[2] * q[4]) * p[3]
    ) +
           (
               (1 - p[2]) * ((q[1] - 1) * f[4] - f[1] * q[4]) +
               ((q[1] - 1) * f[2] - f[1] * q[3]) * p[4] +
               (p[1] - 1) * (f[4] * q[3] - f[2] * q[4])
           ) *
           p[3] *
           q[2] -
           (
               (1 - p[2]) * ((q[1] - 1) * f[3] - (q[2] - 1) * f[1]) +
               ((q[1] - 1) * f[2] - f[1] * q[3]) * p[3] +
               (p[1] - 1) * (f[3] * q[3] - (q[2] - 1) * f[2])
           ) *
           p[4] *
           q[4] -
           (
               (p[1] - 1) * ((q[2] - 1) * f[4] - f[3] * q[4]) +
               ((q[1] - 1) * f[3] - (q[2] - 1) * f[1]) * p[4] -
               ((q[1] - 1) * f[4] - f[1] * q[4]) * p[3]
           ) *
           p[2] *
           q[3])/((p[1] * q[1] - 1) * (
        (q[3] - (q[2] - 1)) * p[4] +
        (p[2] - 1) * ((q[2] - 1) - q[4]) - (q[3] - q[4]) * p[3]
    ) +
           (
               (1 - p[2]) * ((q[1] - 1) - q[4]) +
               ((q[1] - 1) - q[3]) * p[4] +
               (p[1] - 1) * (q[3] - q[4])
           ) *
           p[3] *
           q[2] -
           (
               (1 - p[2]) * ((q[1] - 1) - (q[2] - 1) ) +
               ((q[1] - 1)  -  q[3]) * p[3] +
               (p[1] - 1) * ( q[3] - (q[2] - 1) )
           ) *
           p[4] *
           q[4] -
           (
               (p[1] - 1) * ((q[2] - 1)  -  q[4]) +
               ((q[1] - 1)  - (q[2] - 1) ) * p[4] -
               ((q[1] - 1)  -  q[4]) * p[3]
           ) *
           p[2] *
           q[3])
end



π(p) = π(p, p)

R, S, T, P = RSTP
E(p, q) = (π(p, q) - P) / (π(q, p) - P)
G(p, q) = (R - π(p, q)) / (R - π(q, p))

function extorsion(player, competitors)
    return max(minimum([E(player.strategy, c.strategy) - 1. for c in competitors]), 0.)

end

function generosity(player, competitors)
    return max(minimum([G(player.strategy, c.strategy) - 1. for c in competitors]), 0.)
end

function cooperation(player, competitors)
    return mean([π(player.strategy, c.strategy, [1.0, 1.0, 0.0, 0.0]) for c in competitors])
end

using DrWatson
@quickactivate

using LinearAlgebra
using Test

include(srcdir("IPD-ABM.jl"))

f = rand(4)
p, q = rand(4), rand(4)

@test D(p, q, f)/D(p, q, ones(4)) ≈ dot(v(p, q), f)
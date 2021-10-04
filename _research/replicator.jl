using DrWatson
@quickactivate

using DifferentialEquations
include(srcdir("IPD.jl"))

function replicator!(P, p, t)
    
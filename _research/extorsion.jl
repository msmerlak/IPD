using DrWatson
include(srcdir("memory-one-IPD.jl"))
include(srcdir("constants.jl"))

P = last(ID_PAYOFFS)

EXT = [11/13, 1/2, 7/26, 0.]


q = rand(4)
(π(EXT, q) - P, π(q, EXT) - P)
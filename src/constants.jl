using StaticArrays

RSTP = [3., 0., 5., 1.]

RAND = MVector{4}([.25, .25, .25, .25])

ALLD = MVector{4}([0., 0., 0., 0.])
TFT = MVector{4}([1., 0, 1., 0.])
WSLS = MVector{4}([1., 0., 0., 1.])

PD_EXT = MVector{4}([11/13, 1/2, 7/26, 0.])
SP_GEN = MVector{4}([1., 0.35, .75, .1])

CLASSICAL_STRATEGIES = Vector{MVector{4, Float64}}([TFT, WSLS, PD_EXT, SP_GEN])
using DrWatson
using Dates
using DataFrames
using CSV

include(srcdir("memory-one-IPD.jl"))

X = CSV.read(datadir("mixed", string(today()), "multiplicative=false_n=500_t=10_σ=0.005.csv"), DataFrame)

consensus_strategies = groupby(X, :step)[end]

v
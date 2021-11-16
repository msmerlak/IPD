
using DifferentialEquations, Plots

include(srcdir("memory-one-IPD.jl"))

function ∇π!(dk, x, p, t)
    dk .= ForwardDiff.gradient(π, x)
end

ODEProblem(∇π!, rand(4), (0., 10.), isoutofdomain = (u,p,t)->any(x-> x<0 || x > 1,u)) |> solve |> plot

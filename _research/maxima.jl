using DrWatson
using DifferentialEquations, Plots
using ForwardDiff

include(srcdir("memory-one-IPD.jl"))
include

## fitness maxima
function ∇π!(dk, x, p, t)
    dk .= ForwardDiff.gradient(π, x)
end

ODEProblem(∇π!, rand(4), (0., 10.), isoutofdomain = (u,p,t)->any(x-> x<0 || x > 1,u)) |> solve |> plot



# adaptive dynamics
function f!(dk, x, p, t)
    dk .= ForwardDiff.gradient(y -> π(y, x), x)
end

extorsionary = [11/13, 1/2, 7/26, 0.]

condition(u, t, integrator) = any(x-> x < 0 || x > 1, u) ? 0 : 1
sol = solve(ODEProblem(f!, rand(4), (0., 5.)), callback = ContinuousCallback(condition, terminate!))
plot(sol)

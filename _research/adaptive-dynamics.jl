using DrWatson
using DifferentialEquations, Plots
using ForwardDiff

include(srcdir("memory-one-IPD.jl"))


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


# adaptive dynamics
function f!(dk, x, p, t)
    dk .= ForwardDiff.gradient(y -> π(y, x), x)
end

extorsionary = [11/13, 1/2, 7/26, 0.]

condition(u, t, integrator) = any(x-> x < 0 || x > 1, u) ? 0 : 1
sol = solve(ODEProblem(f!, rand(4), (0., 5.)), callback = ContinuousCallback(condition, terminate!))
plot(sol)


# noisy adaptive dynamics

# function f!(dx, x, p, t)
#     if all(0 .< x .< 1)
#         dx .= ForwardDiff.gradient(y -> π(y, x), x)
#     else
#         @. dx .= -(2x-1)^9
#     end
# end


function f!(dx, x, p, t)
    if all(0 .< x .< 1)
        dx .= ForwardDiff.gradient(y -> π(y, x), x)
    end
    @. dx += -(2x-1)^11
end


function g!(dx, x, p, t)
    @. dx = 5e-2
end

init = rand(4)
T = 20000.
determin = solve(ODEProblem(f!, init, (0., T)))
stoch = solve(SDEProblem(f!, g!, init, (0., T)), maxiters = Int(1e6))
plot(
    plot(determin, ylims = (-.1, 1.1)),
    plot(stoch, ylims = (-0.1, 1.1))
)

plot(stoch.t, map(cooperation, stoch.u), ylims = (-0.1, 1))


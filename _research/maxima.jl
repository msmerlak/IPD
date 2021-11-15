using DrWatson
using IntervalRootFinding

include(srcdir("memory-one-IPD.jl"))

function K((p1, p2, p3, p4))
    return π([p1, p2, p3, p4], [p1, p2, p3, p4])
end

∇K = ∇(K)


using DynamicalSystems

function Q(x, p, t)
    if x ∈ IntervalBox(0.1..0.9, 4)
        return ∇K(x)
    else 
        return zeros(4)
    end
end

gradientflow = ContinuousDynamicalSystem(Q, rand(4), nothing)

B = basins_of_attraction((0:.5:1, 0:.5:1, 0:.5:1, 0:.5:1), gradientflow)

using Optim

lower = zeros(4)
upper = ones(4)
initial_x = rand(4)
inner_optimizer = GradientDescent()
minima = [optimize(x -> -K(x), lower, upper, rand(4), Fminbox(inner_optimizer)).minimizer for _ in 1:100] 
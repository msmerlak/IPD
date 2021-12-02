using DrWatson
@quickactivate

using DifferentialEquations, DiffEqOperators, SparseArrays, LinearAlgebra
using ThreadTools
using Statistics
using Plots


include(srcdir("memory-one-IPD.jl"))

function Δ(X)

	D2 = sparse(CenteredDifference(2, 4, X.step.hi, X.len) * Neumann0BC(X.step.hi))[1]
    id = sparse(I, X.len, X.len)
	
	laplacian = kron(D2, id, id, id) + kron(id, D2, id, id) + kron(id, id, D2, id) + kron(id, id, id, D2)

	return laplacian
end

function replicator!(dP, P, params, t = nothing)
    G = map(x -> sum(P.*x)*params[:dx]^4, params[:Π])
    G .-= mean(G)
    dP .+= P.*G
end 

function mutator!(dP, P, params, t = nothing)
    mul!(dP, params[:Δ], P)
    dP .*= params[:μ]
end

function replicator_mutator!(dP, P, params, t)
    mutator!(dP, P, params, t)
    replicator!(dP, P, params, t)
end

function discretization(X, μ)
    Σ = vec([[p1, p2, p3, p4] for p1 ∈ X, p2 ∈ X, p3 ∈ X, p4 ∈ X])
    Π = tmap(p -> π.(Ref(p), Σ), Σ)
    return Dict(:Π => Π, :Δ => Δ(X), :dx => X.step.hi, :μ => μ)
end

begin
    X = 0.05:.15:.95
    length(X)

    ic = @. exp(-3(X - 0.8)^2)

    pb = ODEProblem(replicator_mutator!, 
        kron(ic, ic, ic, ic), 
        (0., 10.), 
        discretization(X, .1)
        )

    sol = solve(pb, dt = 1.);


    grid = map(x -> reshape(x, (length(X), length(X), length(X), length(X))), sol.u);

end

@gif for t in 1:.1:sol.t[end]
    heatmap(reshape(sol(t), (length(X), length(X), length(X), length(X)))[3, :, :, 3], legend = false, xlabel = "P(C|CC)", title = t)
end



Σ = vec([[p1, p2, p3, p4] for p1 ∈ X, p2 ∈ X, p3 ∈ X, p4 ∈ X])

function mean_strategy(sol, Σ)
    return normalize(sum(sol.*Σ), 1)
end    

plot(transpose([mean_strategy(sol(t), Σ) for t = 0:1:10]))
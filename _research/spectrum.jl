using DrWatson
@quickactivate

using CSV, DataFrames, StatsBase
using Plots
using DSP



# import 10 realizations
F = groupby(CSV.read(datadir("mixed", "one-long-run.csv"), DataFrame), :ensemble)


# focus on one realization
f = F[3].mean_fitness 
plot(f)

ϕ = sign.(f .- 2)

plot(ϕ)

plot(
    F[1].mean_fitness, 
    ylims = (1, 3),
    xlabel = "t",
    ylabel = "mean fitness"
    )

Sf = periodogram(f)
Sϕ = periodogram(ϕ)
plot(Sf.freq[2:end], Sf.power[2:end], xaxis = :log, yaxis = :log)
plot!(Sϕ.freq[2:end], Sϕ.power[2:end], xaxis = :log, yaxis = :log)

plot(Sϕ.freq[2:end], Sϕ.power[2:end], xaxis = :log, yaxis = :log)
plot!(x -> 1e-5/x^2, xlims = (1e-6, 1), yaxis = :log, xaxis = :log)


plot(f, xlims = (1, 50000))
ϕ = round.(f)

x, y = rle(ϕ)
waiting_times = y[x .== 3.]

CDF = ecdf(waiting_times[waiting_times .> 1])

plot(x->1-CDF(x), xlims = (1,maximum(waiting_times)-1), yaxis = :log)

using Distributions
e = fit(Exponential, waiting_times)


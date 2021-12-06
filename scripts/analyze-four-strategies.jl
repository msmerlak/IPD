using DrWatson
@quickactivate

using Arrow, DataFramesMeta, Dates
using StatsPlots, Query


 four_strategies = groupby(DataFrame(
    Arrow.Table(datadir(string(today()),"four-strategies-test.arrow"))), [:init_strategy, :n, :σ])


map(df -> mean(df.mean_cooperation), four_strategies)
combine(four_strategies) do d
    std(d.mean_cooperation)
end

gd = filter(df -> df.n[1] == 10 , four_strategies)

@df combine(gd, [:mean_cooperation => std, :mean_cooperation => mean]) scatter(:σ, :mean_cooperation_mean)
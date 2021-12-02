using DrWatson
@quickactivate

using Arrow, DataFrames, Dates

datadir(string(today()),"four-strategies-n_100.arrow")

four_strategies = Arrow.Table(datadir(string(today()),"four-strategies-n_100.arrow")) |> DataFrame


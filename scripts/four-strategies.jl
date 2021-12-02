using Distributed, Pkg
# @everywhere begin
#     using Pkg
#     Pkg.activate(".")
#     Pkg.instantiate()
# end

PROCESSES = Sys.CPU_THREADS
addprocs(PROCESSES)


@everywhere using DrWatson
@everywhere @quickactivate "IPD"

using Arrow, Dates

### run sims

@everywhere begin
    include(srcdir("IPD.jl"))
    include(srcdir("utils.jl"))
    include(srcdir("constants.jl"))
end


@everywhere function initialize(; n, σ, init_strategy)
    return create_model(
        Dict(:RSTP => RSTP, :n => n, :σ => σ, :init_strategy => init_strategy); 
        compute_metrics = true
        )
end



params = Dict(  :n => 100, 
                :σ => collect(1e-3:1e-3:1e-2), 
                :init_strategy => CLASSICAL_STRATEGIES
            )


adata, _ = myparamscan(params, initialize; 
    agent_step! = player_step!, 
    model_step! = WF_sampling!, 
    n = 500_000, 
    adata = [(:fitness, mean), (:cooperation, mean), (:cooperation, std), (:strategy, mean), (:strategy, std)], 
    parallel = true
    )


try 
    mkdir(datadir(string(today())))
catch
    @warn "today's directory already exists"
end

Arrow.write(datadir(string(today()), "four-strategies-n_100.arrow"), adata)

using Distributed
@everywhere begin
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate()
end

PROCESSES = Sys.CPU_THREADS
addprocs(PROCESSES)



using CSV
CSV.write(datadir("test.csv"), DataFrame(workers = workers()))

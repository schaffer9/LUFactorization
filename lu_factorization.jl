
using LinearAlgebra
using DataFrames
using CSV

include("lu.jl")
include("lu_benchmark_tools.jl")


function run_benchmark_lu_fac()
    # for compilation
    time_solve_lower_tri_sys(10)

    sizes = 100:100:1000

    ts = Float64[]
    factorization_errors = Float64[]
    for n in sizes
        # calculate time for each n
        A, _lu, t = time_lu_fac(n)
        push!(ts, t)
        # calculate errors
        push!(factorization_errors, factorization_error(_lu, A))
    end
    data = DataFrame(size=sizes, time=ts, factorization_error=factorization_errors)

    # save time to csv
    CSV.write("lu_fac.csv", data, writeheaders=true)
end

run_benchmark_lu_fac()

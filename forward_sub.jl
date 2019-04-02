
using LinearAlgebra
using DataFrames
using CSV

include("lu.jl")
include("lu_benchmark_tools.jl")


function run_benchmark_forward_sub()
    # for compilation
    time_solve_lower_tri_sys(10)

    sizes = 100:100:10000

    ts = Float64[]
    forward_errors = Float64[]
    residual_norms = Float64[]
    for n in sizes
        # calculate time for each n
        L, b, x, exact_solution, t = time_solve_lower_tri_sys(n)

        push!(ts, t)
        # calculate errors
        push!(forward_errors, forward_error(L, b, x, exact_solution))
        push!(residual_norms, residual_norm(L, b, x, exact_solution))
    end
    data = DataFrame(size=sizes, time=ts, forward_error=forward_errors, residual_norm=residual_norms)

    # save time to csv
    CSV.write("forward_substitution.csv", data, writeheaders=true)
end

run_benchmark_forward_sub()

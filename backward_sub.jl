
using LinearAlgebra
using DataFrames
using CSV

include("lu.jl")
include("lu_benchmark_tools.jl")


function run_benchmark_backward_sub()
    # for compilation
    time_solve_upper_tri_sys(10)

    sizes = 100:100:10000

    ts = Float64[]
    backward_errors = Float64[]
    residual_norms = Float64[]
    for n in sizes
        # calculate time for each n
        U, b, x, exact_solution, t = time_solve_upper_tri_sys(n)

        push!(ts, t)
        # calculate errors
        push!(backward_errors, forward_error(U, b, x, exact_solution))
        push!(residual_norms, residual_norm(U, b, x, exact_solution))
    end
    data = DataFrame(size=sizes, time=ts, backward_error=backward_errors, residual_norm=residual_norms)

    # save time to csv
    CSV.write("backward_substitution.csv", data, writeheaders=true)
end

run_benchmark_backward_sub()

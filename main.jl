using LUFactorization

function main()
    run_benchmark_lu_fac(100:100:5000, "lu.csv")
    run_benchmark_forward_sub(100:100:5000, "forward_sub.csv")
    run_benchmark_backward_sub(100:100:5000, "backward_sub.csv")
end
main()

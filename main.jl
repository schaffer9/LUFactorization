using LUFactorization
using LinearAlgebra
using CSV

"""
    main(args)

args = [output_path, num_threads, max_size]

output_path is a directory
num_threads is the number of blas threads
max_size is the maximal problem size

```
    julia main.jl . 10 10000
```
"""
function main(args)
    directory_path = isempty(args) ? "." : args[1]
    directory_path = abspath(directory_path)
    mkpath(directory_path)
    
    if length(args) >= 2
        num_blas_thread = parse(Int, args[2])
    else
        num_blas_thread = 1
    end
    BLAS.set_num_threads(num_blas_thread)

    if length(args) >= 3
        max_size = parse(Int, args[3])
    else
        max_size = 3000
    end
    # file_name = "forward_sub.csv"
    # file_path = joinpath(directory_path, file_name)
    # forward_sub_data = run_benchmark_forward_sub(100:100:5000)
    # CSV.write(file_path, forward_sub_data, writeheaders=true)
    #
    # file_name = "backward_sub.csv"
    # file_path = joinpath(directory_path, file_name)
    # backward_sub_data = run_benchmark_backward_sub(100:100:5000)
    # CSV.write(file_path, backward_sub_data, writeheaders=true)

    file_name = "lublas.csv"
    file_path = joinpath(directory_path, file_name)
    lu_with_blas_data = run_benchmark_lu_fac(100:100:min(3000, max_size); blocked=false)
    CSV.write(file_path, lu_with_blas_data, writeheaders=true)

    file_name = "lu_block_size100.csv"
    file_path = joinpath(directory_path, file_name)
    blocked_lu_data100 = run_benchmark_lu_fac(100:100:max_size; blocked=true, block_size=100)
    CSV.write(file_path, blocked_lu_data100, writeheaders=true)

    file_name = "lublocked_different_block_sizes.csv"
    file_path = joinpath(directory_path, file_name)
    blocked_lu_data = run_benchmark_blocked_lu_fac(2000)
    CSV.write(file_path, blocked_lu_data, writeheaders=true)
end

if  !isinteractive()
    main(ARGS)
end

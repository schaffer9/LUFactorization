using DataFrames

factorization_error(factorized::LUFac, A::AbstractMatrix) =
    norm(@view(A[factorized.p, :]) - factorized.L * factorized.U, 1) / norm(@view(A[factorized.p, :]), 1)


function time_lu_fac(n; blocked=true, block_size=100, pivoting=true)
    A = rand(n, n)
    A2 = copy(A)
    sol = @timed lu_factorization!(A2; blocked=blocked, block_size=block_size, pivoting=pivoting)
    return A, sol[1], sol[2]
end


function run_benchmark_lu_fac(sizes; blocked=true, block_size=100, pivoting=true, runs=5)
    # for compilation
    time_lu_fac(100; blocked=blocked, block_size=block_size, pivoting=pivoting)

    ts = Float64[]
    factorization_errors = Float64[]
    for n in sizes
        # calculate time for each n
        t_av = 0
        err_av = 0
        for i = 1:runs
            A, _lu, t = time_lu_fac(n; blocked=blocked, block_size=block_size, pivoting=pivoting)
            t_av += t
            err_av += factorization_error(_lu, A)
        end
        err_av /= runs
        t_av /= runs

        push!(ts, t_av)
        # calculate errors
        push!(factorization_errors, err_av)
    end
    data = DataFrame(size=sizes, time=ts, block_size=fill(block_size, length(sizes)), factorization_error=factorization_errors)

    return data
end

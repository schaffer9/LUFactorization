using DataFrames

factorization_error(factorized::LUFac, A::AbstractMatrix) =
    norm(@view(A[factorized.p, :]) - factorized.L * factorized.U, 1) / norm(@view(A[factorized.p, :]), 1)


forward_error(x, exact_solution) =
    norm(x - exact_solution, 1) / norm(exact_solution, 1)


residual_norm(A, b, x) =
    norm(A * x - b, 1) / norm(b, 1)


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
    forward_err = Float64[]
    res_norm = Float64[]

    for n in sizes
        # calculate time for each n
        t_av = 0.
        fac_err = 0.
        _forward_error = 0.
        _res_norm = 0.
        for i = 1:runs
            A, _lu, t = time_lu_fac(n; blocked=blocked, block_size=block_size, pivoting=pivoting)

            x = ones(n)
            b = A * x
            solution = lu_solve!(_lu, copy(b))
            _forward_error += forward_error(solution, x)
            _res_norm += residual_norm(A, b, solution)
            t_av += t
            fac_err += factorization_error(_lu, A)
        end
        fac_err /= runs
        t_av /= runs
        _forward_error /= runs
        _res_norm /= runs

        push!(ts, t_av)
        push!(forward_err, _forward_error)
        push!(res_norm, _res_norm)
        # calculate errors
        push!(factorization_errors, fac_err)
    end
    data = DataFrame(
        size=sizes,
        time=ts,
        block_size=fill(block_size, length(sizes)),
        factorization_error=factorization_errors,
        forward_error=forward_err,
        residual_norm=res_norm
    )

    return data
end

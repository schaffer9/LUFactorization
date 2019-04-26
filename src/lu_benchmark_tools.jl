using DataFrames

function create_lower_tri_sys(n)
    L = LowerTriangular(rand(n, n))
    L += I * 5
    b = L * ones(n);
    return (L, b)
end

function create_upper_tri_sys(n)
    U = UpperTriangular(rand(n, n))
    U += I * 5
    b = U * ones(n);
    return (U, b)
end

forward_error(A, b, x, exact_solution) =
    norm(x - exact_solution, 1) / norm(exact_solution, 1)

residual_norm(A, b, x, exact_solution) =
    norm(A * x - b, 1) / norm(b, 1)

factorization_error(factorized::LUFac, A::AbstractMatrix) =
    norm(@view(A[factorized.p, :]) - factorized.L * factorized.U, 1) / norm(@view(A[factorized.p, :]), 1)

function factorization_error(A::AbstractMatrix)
    factorized = lu_factorization(A)
    return factorization_error(factorized, A)
end

function time_solve_lower_tri_sys(n)
    L, b = create_lower_tri_sys(n)
    exact_solution = L \ b
    b2 = copy(b)
    sol = @timed forward_substitution!(L, b)
    x = sol[1]
    t = sol[2]
    return (L, b2, x, exact_solution, t)
end

function time_solve_upper_tri_sys(n)
    U, b = create_upper_tri_sys(n)
    exact_solution = U \ b
    b2 = copy(b)
    sol = @timed backward_substitution!(U, b)
    x = sol[1]
    t = sol[2]
    return (U, b2, x, exact_solution, t)
end

function time_lu_fac(n; blocked=true, block_size=100)
    A = rand(n, n)
    A2 = copy(A)
    sol = @timed lu_factorization!(A2; blocked=blocked, block_size=block_size)
    return A, sol[1], sol[2]
end

function run_benchmark_backward_sub(sizes::AbstractArray{<:Int, 1})
    # for compilation
    time_solve_upper_tri_sys(10)

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
    # CSV.write(file_path, data, writeheaders=true)
    return data
end


function run_benchmark_forward_sub(sizes::AbstractArray{<:Int, 1})
    # for compilation
    time_solve_lower_tri_sys(10)

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
    # SV.write(file_path, data, writeheaders=true)
    return data
end

function run_benchmark_lu_fac(sizes::AbstractArray{<:Int, 1}; blocked=true, block_size=100)
    # for compilation
    time_lu_fac(100; blocked=blocked, block_size=block_size)

    ts = Float64[]
    factorization_errors = Float64[]
    for n in sizes
        # calculate time for each n
        t_av = 0
        err_av = 0
        for i = 1:5
            A, _lu, t = time_lu_fac(n; blocked=blocked, block_size=block_size)
            t_av += t
            err_av += factorization_error(_lu, A)
        end
        err_av /= 5
        t_av /= 5

        push!(ts, t_av)
        # calculate errors
        push!(factorization_errors, err_av)
    end
    data = DataFrame(size=sizes, time=ts, factorization_error=factorization_errors)

    # save time to csv
    # CSV.write(file_path, data, writeheaders=true)
    return data
end

function run_benchmark_blocked_lu_fac(s::Integer, block_sizes::AbstractArray{<:Int, 1}=10:10:500)
    # for compilation
    time_lu_fac(100; blocked=true, block_size=100)

    ts = Float64[]
    factorization_errors = Float64[]
    for block_size in block_sizes
        # calculate time for each n
        t_av = 0
        err_av = 0
        for i = 1:20
            A, _lu, t = time_lu_fac(s; blocked=true, block_size=block_size)
            t_av += t
            err_av += factorization_error(_lu, A)
        end
        err_av /= 20
        t_av /= 20
        push!(ts, t_av)
        # calculate errors
        push!(factorization_errors, err_av)
    end
    data = DataFrame(size=fill(s, length(ts)), time=ts, block_size=block_sizes, factorization_error=factorization_errors)

    # save time to csv
    # CSV.write(file_path, data, writeheaders=true)
    return data
end

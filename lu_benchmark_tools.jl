
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


factorization_error(factorized::LU_Fac, A::AbstractMatrix) =
    norm(A - factorized.L * factorized.U, 1) / norm(A, 1)

function factorization_error(A::AbstractMatrix)
    factorized = LU_Fac(A)
    return factorization_error(factorized, A)
end


function time_solve_lower_tri_sys(n)
    L, b = create_lower_tri_sys(n)
    exact_solution = L \ b

    sol = @timed forward_substitution(L, b)
    x = sol[1]
    t = sol[2]
    return (L, b, x, exact_solution, t)
end


function time_solve_upper_tri_sys(n)
    U, b = create_upper_tri_sys(n)
    exact_solution = U \ b

    sol = @timed backward_substitution(U, b)
    x = sol[1]
    t = sol[2]
    return (U, b, x, exact_solution, t)
end


function time_lu_fac(n)
    A = rand(n, n)
    sol = @timed LU_Fac(A)
    return A, sol[1], sol[2]
end

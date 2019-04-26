import LUFactorization:
    forward_substitution,
    forward_substitution!,
    backward_substitution,
    backward_substitution!,
    lu_solve!,
    lu_solve


@testset "forward_substitution" begin
    L = LowerTriangular(rand(1000, 1000)) + I * 10
    b = L * rand(1000)
    x = forward_substitution(L, b)
    @test L * x ≈ b

    L = LowerTriangular(rand(1000, 1000)) + I * 10
    solution = ones(1000)
    b = L * solution
    x = forward_substitution!(L, b)
    @test b ≈ solution

    L = UnitLowerTriangular(rand(1000, 1000))
    solution = rand(1000)
    b = L * solution
    b2 = copy(b)
    x = forward_substitution!(L, b)
    @test L * x ≈ b2

    L = UnitLowerTriangular(rand(1000, 1000))
    solution = ones(1000)
    b = L * solution
    b2 = copy(b)
    x = forward_substitution!(L, b)
    @test L * x ≈ b2
    @test b ≉ b2

end

@testset "backward_substitution" begin
    U = UpperTriangular(rand(1000, 1000)) + I * 10
    b = U * rand(1000)
    x = backward_substitution(U, b)
    @test U * x ≈ b

    U = UpperTriangular(rand(1000, 1000)) + I * 10
    solution = ones(1000)
    b = U * solution
    x = backward_substitution!(U, b)
    @test b ≈ solution

    U = UnitUpperTriangular(rand(1000, 1000))
    solution = rand(1000)
    b = U * solution
    b2 = copy(b)
    x = backward_substitution!(U, b)
    @test U * x ≈ b2

    U = UnitUpperTriangular(rand(1000, 1000))
    solution = ones(1000)
    b = U * solution
    b2 = copy(b)
    x = backward_substitution!(U, b)
    @test U * x ≈ b2
    @test b ≉ b2

end

@testset "Test lu solve" begin
    A = [
        1 2 3 4
        2 3 4 2
        1 1 -1 2
        0 2 1/3 3
    ]
    solution = rand(4, 4)
    b = A * solution
    x = lu_solve(A, b)
    @test A * x ≈ b

    A = rand(1000, 1000)
    A2 = copy(A)
    solution = rand(1000, 1000)
    b = A * solution
    b2 = copy(b)
    x = lu_solve!(A, b)
    @test A2 * b ≈ b2
    @test b ≈ solution

    A = [
        1 2 3 4
        2 3 4 2
        1 1 -1 2
        0 2 1//3 3
    ]
    solution = [1//1, 1, 1, 1]
    b = A * solution
    x = lu_solve(A, b)
    @test A * x ≈ b
    @test x == solution

    A = [
        1 2 3 4
        2 3 4 2
        1 1 -1 2
        0 2 1//3 3
    ]
    solution = [
        1//1 1 1 1
        1 1 1 1
        1 1 1 1
        1 1 1 1
    ]
    b = A * solution
    x = lu_solve(A, b)
    @test A * x ≈ b
    @test x == solution
end

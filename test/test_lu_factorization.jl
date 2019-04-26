import LUFactorization:
    generic_lu_factorization,
    generic_lu_factorization!,
    max_index,
    blocked_lu!,
    blocked_lu


# create LU type
A = rand(100, 100)
F = LUFac(A)

@test F.L == UnitLowerTriangular(A) && F.U == UpperTriangular(A)
@test typeof(F) <: LUFac{<:Real}

A = rand(Int, 5, 5)
F = LUFac(A)

@test F.L == UnitLowerTriangular{Float64}(A) && F.U == UpperTriangular{Float64}(A)
@test typeof(F) <: LUFac{<:Real}


@test F.p == [1,2,3,4,5]

@test propertynames(F) == (:L, :U)
@test propertynames(F, true) == (:lu, :p, :L, :U)


@testset "generic lu factorization" begin
    A = Float64[
        4 3
        6 3
    ]
    L = Float64[
        1 0
        2//3 1
    ]
    U = Float64[
        6 3
        0 1
    ]

    F = generic_lu_factorization(A)
    @test F.L == L && F.U == U

    A = Float64[
        4 3
        6 3
    ]
    L = Float64[
        1 0
        1.5 1
    ]
    U = Float64[
        4 3
        0 -1.5
    ]

    F = generic_lu_factorization(A; pivoting=false)
    @test F.L == L && F.U == U

    A = rand(1000, 1000)
    F2 = generic_lu_factorization(A)
    @test F2.L * F2.U ≈ A[F2.p, :]

    A = rand(Int, 1000, 1000)
    F2 = generic_lu_factorization(A)
    @test F2.L * F2.U ≈ A[F2.p, :]

    # no direct method for Int
    A = rand(Int, 10, 10)
    @test_throws MethodError generic_lu_factorization!(A)

    # matix is overwritten
    A = rand(10, 10)
    F = generic_lu_factorization!(A)
    @test F.L * F.U ≉ A[F.p, :]

    # rational is preserved
    A = [
        1 2 3
        3 1 2
        5 2 1//1
    ]
    F = generic_lu_factorization!(A)
    @test typeof(F) <: LUFac{<:Rational}
end

@testset "find max elem in column" begin
    A = [
        1 2 3
        3 1 2
        5 2 1
    ]
    @test max_index(A, 1) == 3
    @test max_index(A, 2) == 3
    @test max_index(A, 3) == 3
end

@testset "Test blocked version" begin
    A = [
        4 3 4 5 2 1
        3 4 8 1 2 7
        7 8 9 2 1 3
        7 1 2 3 4 5
        1 9 7 3 3 1
        6 5 6 1 2 9.
    ]
    F1 = generic_lu_factorization(A; pivoting=false)
    F2 = blocked_lu(A, 3; pivoting=false)
    @test F2.L * F2.U ≈ A

    A = rand(1000, 1000)

    F = blocked_lu(A, 100)
    @test F.L * F.U ≈ A[F.p, :]
end


@testset "Test lu factorization" begin
    A = rand(1000, 1000)
    F = lu_factorization(A)
    @test F.L * F.U ≈ A[F.p, :]

    A = rand(1000, 1000)
    F = lu_factorization(A; blocked=false)
    @test F.L * F.U ≈ A[F.p, :]
end

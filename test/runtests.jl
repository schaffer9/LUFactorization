using LUFactorization
using LinearAlgebra
using Test

@testset "Test LUFac" begin include("test_lu_fac.jl") end
@testset "Test basic lu factorization" begin include("test_lu_factorization.jl") end
@testset "Test substitution" begin include("test_substitution.jl") end

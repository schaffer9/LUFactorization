using LUFactorization
using LinearAlgebra
using Test

@testset "Test basic lu factorization" begin include("test_lu_factorization.jl") end
@testset "Test substitution" begin include("test_substitution.jl") end

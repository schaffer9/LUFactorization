
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

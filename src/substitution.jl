import LinearAlgebra: ldiv!
import Base: \

function forward_substitution!(L::LowerTriangular, b::AbstractArray)
    b[1] = b[1] / L[1,1]
    for k = 2:length(b)
        for j = 1:k-1
            b[k] -= L[k,j] * b[j]
        end
        b[k] /= L[k,k]
    end
    return b
end


function backward_substitution!(U::UpperTriangular, b::AbstractArray)
    l = length(b)
    b[end] = b[end] / U[end,end]
    for k = (l - 1):-1:1
        for j = k+1:l
            b[k] -= U[k,j] * b[j]
        end
        b[k] /= U[k,k]
    end
    return b
end


function forward_substitution!(L::UnitLowerTriangular, b::AbstractArray)
    for k = 2:length(b)
        for j = 1:k-1
            b[k] -= L[k,j] * b[j]
        end
    end
    return b
end


function backward_substitution!(U::UnitUpperTriangular, b::AbstractArray)
    l = length(b)
    for k = (l - 1):-1:1
        for j = k+1:l
            b[k] -= U[k,j] * b[j]
        end
    end
    return b
end

forward_substitution(L::AbstractMatrix, b::AbstractArray) = forward_substitution!(L, copy(b))
backward_substitution(U::AbstractMatrix, b::AbstractArray) = backward_substitution!(U, copy(b))

function _solve!(F::LUFac{T}, b::AbstractArray{T}) where T<:AbstractFloat
    BLAS.trsv!('L', 'N', 'U', F.lu, bp)
    BLAS.trsv!('U', 'N', 'N', F.lu, bp)
end


function _solve!(F::LUFac{T}, B::AbstractMatrix{T}) where T<:AbstractFloat
    BLAS.trsm!('L', 'L', 'N', 'U', 1.0, F.lu, B)
    BLAS.trsm!('L', 'U', 'N', 'N', 1.0, F.lu, B)
end

function _solve!(F::LUFac{<:Real}, b::AbstractArray{<:Real})
    forward_substitution!(UnitLowerTriangular(F.lu), b)
    backward_substitution!(UpperTriangular(F.lu), b)
end

function _solve!(A::LUFac, B::AbstractMatrix{<:Real})
    mapslices(B; dims=1) do b _solve!(A, b) end
end

function swap_pivot_rows!(b::AbstractArray, piv::AbstractArray{<:Integer, 1})
    bp = b[piv, :]
    copyto!(b, bp)
end

function lu_solve!(lu::LUFac{T}, b::AbstractArray{T}) where T<:Real
    swap_pivot_rows!(b, lu.p)
    _solve!(lu, b)
end

function lu_solve!(A::AbstractMatrix{T}, b::AbstractArray{T}) where T<:Real
    lu = lu_factorization!(A)
    lu_solve!(lu, b)
end

lu_solve(A::Union{AbstractMatrix, LUFac}, b::AbstractArray) = lu_solve!(copy(A), copy(b))

\(A::LUFac, b::AbstractArray) = lu_solve(A, b)
ldiv!(A::LUFac, b::AbstractArray) = lu_solve!(A, b);


function forward_substitution(L, b)
    x = similar(b)

    x[1] = b[1] / L[1,1]
    @inbounds begin
        for k = 2:length(b)
            x[k] = b[k]

            for j = 1:k-1
                x[k] -= L[k,j] * x[j]
            end

            x[k] /= L[k,k]
        end
    end
    return x
end


function backward_substitution(L, b)
    l = length(b)
    x = similar(b)

    x[end] = b[end] / L[end,end]
    @inbounds begin
        for k = (l - 1):-1:1
            x[k] = b[k]

            for j = k+1:l
                x[k] -= L[k,j] * x[j]
            end
            x[k] /= L[k,k]
        end
    end
    return x
end

function lu_factorization(A::AbstractMatrix)
    n, m = size(A)
    @assert n == m
    L = Matrix{Float64}(I, n, n)
    U = copy(A)
    @inbounds begin
        for k = 1:n-1
                # L[j,k] = U[j,k] / U[k,k]
            Ukkinv = inv(U[k,k])
            for j = k+1:n
                L[j,k] = U[j,k] * Ukkinv  # U[k,k]
                for l = k:n
                    U[j,l] = U[j,l] - L[j,k] * U[k,l]
                end
            end
        end
    end

    return (L=L, U=U)
end

struct LU_Fac
    L::AbstractMatrix
    U::AbstractMatrix
    function LU_Fac(A::AbstractMatrix)
        lu = lu_factorization(A)
        new(lu.L, lu.U)
    end
end


function lu_solve(A::LU_Fac, b)
    y = forward_substitution(A.L, b)
    x = backward_substitution(A.U, y)
end

function lu_solve(A::AbstractArray, b)
    lu = LU_Fac(A)
    y = forward_substitution(lu.L, b)
    x = backward_substitution(lu.U, y)
end

import Base: \
\(A::LU_Fac, b::AbstractArray) = lu_solve(A, b);

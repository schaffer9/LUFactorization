
@inline max_index(A::AbstractMatrix, col::Integer) = findmax(@view(A[col:end, col]))[2] + col - 1


@inline function swap_rows!(A::AbstractArray{<:Real}, i, j)
    A[i], A[j] = A[j], A[i]
    return A
end


@inline function swap_rows!(A::AbstractMatrix{<:Real}, i, j)
    for l = 1:size(A, 2)
        A[i, l], A[j, l] = A[j, l], A[i, l]
    end
    return A
end


@inline function _factorize!(A::AbstractMatrix{<:Real}, k::Integer, n::Integer, m::Integer)
    Akkinv = inv(A[k,k])
    @inbounds begin
        for j = k+1:n
            A[j,k] *= Akkinv
        end
        for l = k+1:m
            for j = k+1:n
                A[j,l] -= A[j,k] * A[k,l]
            end
        end
    end
end

@inline function _factorize!(A::AbstractMatrix{<:AbstractFloat}, k::Integer, n::Integer, m::Integer)
    Akkinv = inv(A[k,k])
    l = @view(A[k+1:n,k])
    BLAS.scal!(length(l), Akkinv, l, 1)
    @views BLAS.ger!(-1.0, A[k+1:n,k], A[k,k+1:m], A[k+1:n,k+1:m])
end


function _lu_factorization!(
    A::AbstractMatrix,
    piv::AbstractArray{<:Integer},
    i::OrdinalRange,
    j::OrdinalRange;
    pivoting=true
)
    Aij = @view A[i, j]
    n, m = size(Aij)
    @inbounds begin
        for k = 1:m
            if pivoting
                kp = k + i.start - 1
                pivot = max_index(A, kp)
                if kp != pivot
                    swap_rows!(piv, kp, pivot)
                    swap_rows!(A, kp, pivot)
                end
            end
            _factorize!(Aij, k, n, m)
        end
    end
    A
end


function generic_lu_factorization!(A::AbstractMatrix{T}; pivoting=true) where T<:Real
    if T <: Integer
        throw(MethodError("Cannot factorize Integer matrix in place!"))
    end
    n, m = size(A)
    if m > n
        error("Dimension mismatch! Cannot factorize matrix with more columns than rows.")
    end
    piv = collect(1:n)

    _lu_factorization!(A, piv, 1:n, 1:m; pivoting=pivoting)

    return LUFac{T}(
        A,
        piv
    )
end

generic_lu_factorization(A; pivoting=true) = generic_lu_factorization!(copy(A); pivoting=pivoting)
generic_lu_factorization(A::AbstractMatrix{T}; pivoting=true) where T<:Integer =
    generic_lu_factorization(convert(Matrix{Float64}, A); pivoting=pivoting)


_last_block(n, b) = n%b == 0 ? n - b +1 : n - (n%b) + 1

function blocked_lu!(A::AbstractMatrix{T}, b::Integer; pivoting=true) where T <: AbstractFloat
    n = min(size(A)...)
    piv = collect(1:n)
    last_block = _last_block(n, b)
    @inbounds begin
        for i = 1:b:last_block-b
            k = i:i+b-1
            l = i+b:n
            _lu_factorization!(A, piv, i:n, k; pivoting=pivoting)
            Lkk = @view A[k, k]
            Ukl = @view A[k, l]
            Llk = @view A[l, k]
            All = @view A[l,  l]
            # Aji = Lkk \ Aji
            BLAS.trsm!('L', 'L', 'N', 'U', 1.0, Lkk, Ukl)
            # All = All - Llk*Ukl
            BLAS.gemm!('N', 'N', -1.0, Llk, Ukl, 1.0, All)
        end
    end
    _lu_factorization!(A, piv, last_block:n, last_block:n; pivoting=pivoting)
    return LUFac{T}(
        A,
        piv
    )
end

blocked_lu(A::AbstractMatrix, b::Integer; pivoting=true) = blocked_lu!(copy(A), b; pivoting=pivoting)

function lu_factorization!(A::AbstractMatrix{T}; pivoting::Bool=true, blocked::Bool=true, block_size::Integer=100) where T
    if blocked && T<:AbstractFloat
        return blocked_lu!(A, block_size; pivoting=pivoting)
    else
        return generic_lu_factorization!(A; pivoting=pivoting)
    end
end

lu_factorization(A::AbstractMatrix; pivoting::Bool=true, blocked::Bool=true, block_size::Integer=100) =
    lu_factorization!(copy(A); pivoting=pivoting, blocked=blocked, block_size=block_size)

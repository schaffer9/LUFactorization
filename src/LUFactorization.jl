module LUFactorization

using LinearAlgebra
import Base: getproperty, propertynames, show

struct LUFac{T<:Real}
    lu::Matrix{T}
    p::Array{Int, 1}

    LUFac{T}(lu::AbstractArray) where T = LUFac{T}(lu, collect(1:size(lu, 1)))
    LUFac{T}(lu::AbstractArray, p::AbstractArray{<:Integer, 1}) where T<:Real =
        new{T}(convert(Matrix{T}, lu), p)

    LUFac(lu::AbstractArray)= LUFac(lu, collect(1:size(lu, 1)))
    LUFac(lu::AbstractArray, p::AbstractArray{<:Integer, 1}) = LUFac{Float64}(lu, p)

end

export LUFac

function getproperty(F::LUFac, d::Symbol)

    if d === :L
        return UnitLowerTriangular(F.lu)
    elseif d === :U
        return UpperTriangular(F.lu)
    else
        getfield(F, d)
    end
end

function propertynames(F::LUFac, private::Bool=false)
    properties = (:L, :U)
    if private
        return (fieldnames(typeof(F))..., properties...)
    else
        return properties
    end
end

function show(io::IO, mime::MIME{Symbol("text/plain")}, F::LUFac)
    print(io, "L = ")
    show(io, mime, F.L)
    print(io, "\n\nU = ")
    show(io, mime, F.U)
end

include("lu_factorization.jl")
export lu_factorization, lu_factorization!


include("substitution.jl")
export lu_solve, lu_solve!

include("lu_benchmark_tools.jl")
export run_benchmark_lu_fac
end # module

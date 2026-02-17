module ArrayInterfaceDynamicQuantitiesExt

using ArrayInterface
using DynamicQuantities
using LinearAlgebra

# DynamicQuantities quantities have runtime-stored dimensions, so identity elements
# like `zero(::Type{<:Quantity})` / `oneunit(::Type{<:Quantity})` are intentionally
# undefined. ArrayInterface's instance constructors should therefore avoid calling
# `zeros(T, ...)` on Quantity element types.

# Avoid method ambiguities with ArrayInterface's potential `T<:Number` specializations
# by defining these methods per-abstract quantity type, rather than using
# `T<:UnionAbstractQuantity` (a `Union` bound).
for (Qabs, _, _) in DynamicQuantities.ABSTRACT_QUANTITY_TYPES
    @eval begin
        function ArrayInterface.qr_instance(
            A::Matrix{T},
            pivot = ArrayInterface.DEFAULT_CHOLESKY_PIVOT,
        ) where {T<:$Qabs}
            if pivot === ArrayInterface.DEFAULT_CHOLESKY_PIVOT
                return LinearAlgebra.QRCompactWY(Matrix{T}(undef, 0, 0), Matrix{T}(undef, 0, 0))
            else
                return LinearAlgebra.QRPivoted(
                    Matrix{T}(undef, 0, 0),
                    Vector{T}(undef, 0),
                    Vector{Int}(undef, 0),
                )
            end
        end

        function ArrayInterface.svd_instance(A::Matrix{T}) where {T<:$Qabs}
            return LinearAlgebra.SVD(
                Matrix{T}(undef, 0, 0),
                Vector{real(T)}(undef, 0),
                Matrix{T}(undef, 0, 0),
            )
        end

        # Provide a cheap LU "instance" for Quantity matrices without calling `lu(A)`.
        # This is used for cache allocation; actual factorization should be done elsewhere.
        function ArrayInterface.lu_instance(A::Matrix{T}) where {T<:$Qabs}
            isempty(A) && return ArrayInterface.lu_instance(Matrix{Float64}(undef, 0, 0))
            noUnitT = typeof(zero(DynamicQuantities.ustrip(first(A))))
            luT = LinearAlgebra.lutype(noUnitT)
            ipiv = Vector{LinearAlgebra.BlasInt}(undef, 0)
            info = zero(LinearAlgebra.BlasInt)
            return LinearAlgebra.LU{luT}(similar(A, 0, 0), ipiv, info)
        end
    end
end

end

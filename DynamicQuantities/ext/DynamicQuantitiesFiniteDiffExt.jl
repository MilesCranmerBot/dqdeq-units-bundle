module DynamicQuantitiesFiniteDiffExt

using DynamicQuantities
import FiniteDiff

# FiniteDiff’s default step-size logic uses `real(T)` for the scalar type `T`.
# For DynamicQuantities quantities, `real(::Type{<:Quantity})` falls back to
# `zero(T)` which is intentionally undefined (dims are runtime).
# We can delegate to the underlying primitive numeric type instead.
FiniteDiff.default_relstep(v::Val, ::Type{<:UnionAbstractQuantity{T}}) where {T} =
    FiniteDiff.default_relstep(v, T)

end

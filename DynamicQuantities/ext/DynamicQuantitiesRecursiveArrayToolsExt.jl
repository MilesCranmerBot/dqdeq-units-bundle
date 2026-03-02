module DynamicQuantitiesRecursiveArrayToolsExt

using DynamicQuantities
using RecursiveArrayTools

# Ensure RecursiveArrayTools' unitless type extraction strips
# DynamicQuantities wrappers instead of stopping at the quantity
# type itself. This keeps existing `recursive_unitless_*` behavior
# compatible with value-space assumptions in DifferentialEquations paths.
for (Qabs, _, _) in DynamicQuantities.ABSTRACT_QUANTITY_TYPES
    @eval begin
        RecursiveArrayTools.recursive_unitless_bottom_eltype(::Type{<:$Qabs{T,D}}) where {T,D} =
            RecursiveArrayTools.recursive_unitless_bottom_eltype(T)
        RecursiveArrayTools.recursive_unitless_eltype(::Type{<:$Qabs{T,D}}) where {T,D} =
            RecursiveArrayTools.recursive_unitless_eltype(T)
    end
end

end

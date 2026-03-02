module DynamicQuantitiesSciMLBaseExt

using DynamicQuantities
import SciMLBase: unitfulvalue, value

# SciMLBase defines `value`/`unitfulvalue` hooks for unitful scalar types.
# Define these for DynamicQuantities quantities here (owner of the quantity types)
# to avoid type piracy in downstream SciML packages.
for (_Q, _, _) in DynamicQuantities.ABSTRACT_QUANTITY_TYPES
    @eval begin
        value(::Type{<:$_Q{T}}) where {T} = T
        value(x::$_Q) = ustrip(x)

        unitfulvalue(::Type{T}) where {T<:$_Q} = T
        unitfulvalue(x::$_Q) = x
    end
end

end

function pre_basis(a::AbstractArray{T}) where {T}
    # Prefer value-based identities so scalar types without `zero(::Type{T})`
    # (e.g. runtime-unit quantities) can still participate.
    #
    # Arrays may (in principle) contain mixed "units" (e.g. runtime-unit quantities
    # or heterogeneous containers), so construct the zero array elementwise.
    out = similar(a)
    @inbounds for i in eachindex(a)
        out[i] = zero(a[i])
    end
    return out
end

function post_basis(b::AbstractArray, a::AbstractArray)
    if ismutable_array(a)
        return b
    else
        return map(+, zero(a), b)
    end
end

"""
    basis(a::AbstractArray, i)

Construct the `i`-th standard basis array in the vector space of `a`.
"""
function basis(a::AbstractArray, i)
    b = pre_basis(a)
    b[i] = oneunit(a[i])
    return post_basis(b, a)
end

# compatible with zero-length vectors
function basis(a::AbstractArray)
    b = pre_basis(a)
    return post_basis(b, a)
end

"""
    multibasis(a::AbstractArray, inds)

Construct the sum of the `i`-th standard basis arrays in the vector space of `a` for all `i ∈ inds`.
"""
function multibasis(a::AbstractArray, inds)
    b = pre_basis(a)
    isempty(b) && return post_basis(b, a)
    for i in inds
        b[i] = oneunit(a[i])
    end
    return post_basis(b, a)
end

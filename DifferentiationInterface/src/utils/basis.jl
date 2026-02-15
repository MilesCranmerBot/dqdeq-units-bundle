function pre_basis(a::AbstractArray{T}) where {T}
    # Prefer value-based identities so scalar types without `zero(::Type{T})`
    # (e.g. runtime-unit quantities) can still participate.
    isempty(a) && return similar(a)
    return fill!(similar(a), zero(first(a)))
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
    b[i] = isempty(b) ? oneunit(eltype(b)) : oneunit(first(b))
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
    oneb = isempty(b) ? oneunit(eltype(b)) : oneunit(first(b))
    for i in inds
        b[i] = oneb
    end
    return post_basis(b, a)
end

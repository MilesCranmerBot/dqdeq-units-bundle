module DiffEqBaseDynamicQuantitiesExt

using DiffEqBase
import SciMLBase: unitfulvalue, value
using DynamicQuantities
using LinearAlgebra
import DiffEqBase: default_factorize

# Support adaptive errors should be errorless for exponentiation
value(::Type{<:UnionAbstractQuantity{T}}) where {T} = T
value(x::UnionAbstractQuantity) = ustrip(x)

unitfulvalue(::Type{T}) where {T <: UnionAbstractQuantity} = T
unitfulvalue(x::UnionAbstractQuantity) = x

@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{<:UnionAbstractQuantity, N},
        t,
    ) where {N}
    return sqrt(
        sum(
            x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
            zip((value(x) for x in u), Iterators.repeated(t)),
        ) / length(u),
    )
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::Array{<:UnionAbstractQuantity, N},
        t,
    ) where {N}
    return sqrt(
        sum(
            x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
            zip((value(x) for x in u), Iterators.repeated(t)),
        ) / length(u),
    )
end
@inline DiffEqBase.ODE_DEFAULT_NORM(u::UnionAbstractQuantity, t) = abs(value(u))
@inline function DiffEqBase.UNITLESS_ABS2(x::UnionAbstractQuantity)
    return real(abs2(ustrip(x)))
end

DiffEqBase._rate_prototype(u, t::UnionAbstractQuantity, onet) = u / oneunit(t)
DiffEqBase.timedepentdtmin(t::UnionAbstractQuantity, dtmin) =
    abs(value(dtmin / oneunit(t)) * oneunit(t))

# Rosenbrock/SDIRK solvers form W/J matrices with Quantity eltype. Factorize/solve in
# value-space (Float64), but return solutions with the RHS units.
struct DQUnitlessLU{F, UT}
    F::F
    ut::UT
end

@inline function _infer_ut(A::AbstractMatrix{<:UnionAbstractQuantity})
    @inbounds for a in A
        va = value(a)
        if !iszero(va)
            return oneunit(inv(a))
        end
    end
    return oneunit(1.0)
end

function default_factorize(A::AbstractMatrix{<:UnionAbstractQuantity})
    isempty(A) && return DQUnitlessLU(
        lu(Matrix{Float64}(undef, 0, 0); check = false),
        oneunit(1.0),
    )
    ut = _infer_ut(A)
    return DQUnitlessLU(lu(value.(A); check = false), ut)
end

function LinearAlgebra.ldiv!(
        x::AbstractVector{<:UnionAbstractQuantity},
        W::DQUnitlessLU,
        b::AbstractVector{<:UnionAbstractQuantity},
    )
    vb = value.(b)
    vx = similar(vb)
    LinearAlgebra.ldiv!(vx, W.F, vb)
    @inbounds for i in eachindex(x)
        x[i] = vx[i] * (oneunit(b[i]) * W.ut)
    end
    return x
end

function Base.:(\)(W::DQUnitlessLU, b::AbstractVector{<:UnionAbstractQuantity})
    vb = value.(b)
    vx = W.F \ vb
    out = similar(b)
    @inbounds for i in eachindex(out)
        out[i] = vx[i] * (oneunit(b[i]) * W.ut)
    end
    return out
end

end

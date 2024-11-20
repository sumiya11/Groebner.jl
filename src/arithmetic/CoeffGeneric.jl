# This file is a part of Groebner.jl. License is GNU GPL v2.

struct CoeffGeneric{T}
    data::T
end

generic_coeff_data(x) = x
generic_coeff_data(x::CoeffGeneric) = x.data

CoeffGeneric(x::CoeffGeneric{T}) where {T} = x
CoeffGeneric{T}(x::CoeffGeneric{T}) where {T} = x

Base.:+(x::CoeffGeneric, y::CoeffGeneric) = CoeffGeneric(x.data + y.data)
Base.:-(x::CoeffGeneric, y::CoeffGeneric) = CoeffGeneric(x.data - y.data)
Base.:*(x::CoeffGeneric, y::CoeffGeneric) = CoeffGeneric(x.data * y.data)
Base.inv(x::CoeffGeneric) = CoeffGeneric(inv(x.data))
Base.:-(x::CoeffGeneric) = zero(x) - x

Base.one(x::CoeffGeneric{T}) where {T} = CoeffGeneric{T}(one(x.data))
Base.zero(x::CoeffGeneric{T}) where {T} = CoeffGeneric{T}(zero(x.data))

Base.isone(x::CoeffGeneric) = isone(x.data)
Base.iszero(x::CoeffGeneric) = iszero(x.data)

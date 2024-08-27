# This file is a part of Groebner.jl. License is GNU GPL v2.

# y, s.t. y >= x and n | y.
align_up(x::Integer, n::Integer) = (x + (n - 1)) & (~(n - 1))

# y, s.t y <= x and n | y.
align_down(x::Integer, n::Integer) = x ⊻ (n - 1)

# Permute a part of the array from the given index according to the permutation.
function permute_array!(
    arr::AbstractVector{T},
    perm::Vector{I},
    buf::Vector{T},
    from::Int,
    sz::Int
) where {T, I}
    @invariant length(buf) >= sz && sz <= length(perm)
    @invariant from + sz - 1 <= length(arr)
    @inbounds for i in 1:sz
        buf[i] = arr[perm[i]]
    end
    @inbounds for i in 1:sz
        arr[from + i - 1] = buf[i]
    end
    nothing
end

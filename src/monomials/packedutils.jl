# Utilities for packed integers

# Functions here work on integers packed into a single integer of a wider type.

# How many integers of type B can be stored in an integer of type T
function packed_elperchunk(T, B)
    epc = div(sizeof(T), sizeof(B))
    @invariant epc * sizeof(B) == sizeof(T)
    epc
end

# the size of total degree.
# By default, does not differ from the size of other elements
packed_degsize(T, B, n) = sizeof(B)
packed_nchunks(T, B, n) = div((n - 1) * sizeof(B) + packed_degsize(T, B, n), sizeof(T)) + 1

@generated function packed_unpack!(
    b::AbstractVector{MH},
    a::T,
    ::Type{B},
    I::Integer
) where {T, MH, B}
    ts, bs = sizeof(T), sizeof(B)
    @assert bs * div(ts, bs) == ts
    epc = div(ts, bs)
    shift = bs * 8
    ans = :()
    for i in 1:epc
        ans = :($ans;
        ($i + I > $epc) && return $b;
        @inbounds b[$i] = mod(a, $B);
        a = a >> $shift)
    end
    :($ans; return $b)
end

@generated function packed_dot_product(
    a::T,
    b::AbstractVector{MH},
    ::Type{B},
    I::Integer
) where {T, MH, B}
    ts, bs = sizeof(T), sizeof(B)
    @assert bs * div(ts, bs) == ts
    epc = div(ts, bs)
    shift = bs * 8
    x = :x
    ans = :(x = $MH(0))
    for i in 1:epc
        ans = :($ans;
        ($i + I > $epc) && return $x;
        @inbounds iz = $MH(mod(a, $B)) * b[$i];
        x = x + iz;
        a = a >> $shift)
    end
    :($ans; return $x)
end

# Keep the number of Val(...) parametric types small, otherwise these functions
# will be dispatched in runtime
@generated function packed_is_zero_dot_product(
    a::T,
    b::T,
    ::Type{B},
    ::Val{I}
) where {T, B, I}
    ts, bs = sizeof(T), sizeof(B)
    @assert bs * div(ts, bs) == ts
    epc = div(ts, bs)
    shift = bs * 8
    x = :x
    ans = :(x = true)
    for i in 1:(epc - I)
        ans = :($ans;
        iz = iszero(mod(a, $B)) || iszero(mod(b, $B));
        x = $x && iz;
        a = a >> $shift;
        b = b >> $shift)
    end
    :($ans; return $x)
end

# a >= b
@generated function packed_ge(a::T, b::T, ::Type{B}, ::Val{I}) where {T, B, I}
    ts, bs = sizeof(T), sizeof(B)
    @assert bs * div(ts, bs) == ts
    epc = div(ts, bs)
    shift = bs * 8
    x = :x
    ans = :(a < b && return false; x = true)
    for i in 1:(epc - I)
        ans = :($ans;
        iz = mod(a, $B) >= mod(b, $B);
        x = $x && iz;
        a = a >> $shift;
        b = b >> $shift)
    end
    :($ans; return $x)
end

@generated function packed_max(a::T, b::T, ::Type{B}, ::Val{I}) where {T, B, I}
    ts, bs = sizeof(T), sizeof(B)
    @assert bs * div(ts, bs) == ts
    epc = div(ts, bs)
    shift = bs * 8
    x = :x
    si = :si
    ans = :(x = zero($T); si = zero($T))
    for i in 1:(epc - I)
        ans = :($ans;
        iz = $T(max(mod(a, $B), mod(b, $B)));
        x = x | (iz << ($shift * ($i - 1)));
        si = si + iz;
        a = a >> $shift;
        b = b >> $shift)
    end
    :($ans; return $x, $si)
end
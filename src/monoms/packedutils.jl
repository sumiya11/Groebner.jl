# Utilities for packed monomial representation

# Threshold for overflow error to trigger.
# > If entries in exponent vector are of type B,
# > then typemax(B)/2 is enough.
# That is, the maximal element of the exponent vector is the total degree.
# The total degree can grow with
# monom-monom multiplications (no more than twice)
# and monom-monom lcm (no more than twice). 
_overflow_threshold(B) = div(typemax(B), 2)
@noinline _overflow_error(c, B) = throw(RecoverableException("Overflow is probable with $c modulo $B."))

function _overflow_check(e::Integer, B)
    e >= _overflow_threshold(B) && _overflow_error(e, B)
    true
end

# How many integers of type B can be stored
# in an integer of type T
function elperchunk(T, B) 
    epc = div(sizeof(T), sizeof(B))
    @assert epc*sizeof(B) == sizeof(T)
    epc
end

# the size of total degree.
# By default, does not differ from the size of other elements
degsize(T, B, n) = sizeof(B)
nchunks(T, B, n) = div((n-1)*sizeof(B) + degsize(T, B, n), sizeof(T)) + 1

# For an integer a::T that packs several integers of type B
# unpacks and writes the stored integers into vector `b`.
#
# keep the number of Val(...) parametric types small,
# otherwise these functions will be dispatched in runtime
@generated function packedunpack!(b::AbstractVector{MH}, a::T, ::Type{B}, I::Integer) where {T, MH, B}
    ts, bs = sizeof(T), sizeof(B)
    @assert bs*div(ts, bs) == ts
    epc = div(ts, bs)
    shift = bs*8
    ans = :()
    for i in 1:epc
        ans = :(    
            $ans;
            ($i + I > $epc) && return $b;
            @inbounds b[$i] = mod(a, $B);
            a = a >> $shift;
        )
    end
    :($ans; return $b)
end

# For an integer a::T that packs several integers of type B
# computes the dot product of a and vector `b`
@generated function packeddot(a::T, b::AbstractVector{MH}, ::Type{B}, I::Integer) where {T, MH, B}
    ts, bs = sizeof(T), sizeof(B)
    @assert bs*div(ts, bs) == ts
    epc = div(ts, bs)
    shift = bs*8
    x = :x
    ans = :(x = $MH(0);)
    for i in 1:epc
        ans = :(
            $ans;
            ($i + I > $epc) && return $x;
            @inbounds iz = $MH(mod(a, $B)) * b[$i];
            x = x + iz;
            a = a >> $shift;
        )
    end
    :($ans; return $x)
end

# For integers a::T and b::T that pack several integers of type B
# checks that dot product of a and b is zero
@generated function packedorth(a::T, b::T, ::Type{B}, ::Val{I}) where {T, B, I}
    ts, bs = sizeof(T), sizeof(B)
    @assert bs*div(ts, bs) == ts
    epc = div(ts, bs)
    shift = bs*8
    x = :x
    ans = :(x = true;)
    for i in 1:epc-I
        ans = :(
            $ans; 
            iz = iszero(mod(a, $B)) || iszero(mod(b, $B));
            x = $x && iz;
            a = a >> $shift;
            b = b >> $shift;
        )
    end
    :($ans; return $x)
end

# For integers a::T and b::T that pack several integers of type B
# checks that a >= b lexicographically
@generated function packedge(a::T, b::T, ::Type{B}, ::Val{I}) where {T, B, I}
    ts, bs = sizeof(T), sizeof(B)
    @assert bs*div(ts, bs) == ts
    epc = div(ts, bs)
    shift = bs*8
    x = :x
    ans = :(a < b && return false; x = true;)
    for i in 1:epc-I
        ans = :(
            $ans; 
            iz = mod(a, $B) >= mod(b, $B);
            x = $x && iz;
            a = a >> $shift;
            b = b >> $shift;
        )
    end
    :($ans; return $x)
end

# For integers a::T and b::T that pack several integers of type B
# computes and returns
# max.(a, b), sum(max.(a, b))
@generated function packedmax(a::T, b::T, ::Type{B}, ::Val{I}) where {T, B, I}
    ts, bs = sizeof(T), sizeof(B)
    @assert bs*div(ts, bs) == ts
    epc = div(ts, bs)
    shift = bs*8
    x = :x
    si = :si
    ans = :(x = zero($T); si = zero($T);)
    for i in 1:epc-I
        ans = :(
            $ans; 
            iz = $T(max(mod(a, $B), mod(b, $B)));
            x = x | (iz << ($shift*($i-1))); si = si + iz;
            a = a >> $shift;
            b = b >> $shift;
        )
    end
    :($ans; return $x, $si)
end

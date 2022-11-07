
_overflow_threshold(B) = div(typemax(B), 2)
@noinline _overflow_error(c, B) = throw(OverflowError("Overflow is probable with $c modulo $B."))

function _overflow_check(e::Integer, B)
    e >= _overflow_threshold(B) && _overflow_error(e, B)
    true
end

function elperchunk(T, B) 
    epc = div(sizeof(T), sizeof(B))
    @assert epc*sizeof(B) == sizeof(T)
    epc
end

degsize(T, B, n) = sizeof(B)
nchunks(T, B, n) = div((n-1)*sizeof(B) + degsize(T, B, n), sizeof(T)) + 1

# keep the number of parametric types small,
# otherwise these functions could be dispatched in runtime
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
            b[$i] = mod(a, $B);
            a = a >> $shift;
        )
    end
    :($ans; return $b)
end

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

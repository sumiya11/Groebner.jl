# Reduction of the lower rows B at a time: the accumulators of the B rows are
# interleaved, A[(c - 1) * B + l] holds column c of lane l, in [0, p^2).

# a[0..B) -= m[0..B) * x, kept in [0, p^2): one vector update, products 32 x 32 -> 64
@generated function _blk_vupd!(
    pa::Ptr{Int64},
    m::NTuple{B, VecElement{UInt32}},
    x::UInt32,
    p2::Int64
) where {B}
    ir = """
    %a = load <$B x i64>, ptr %0, align 8
    %mz = zext <$B x i32> %1 to <$B x i64>
    %xz = zext i32 %2 to i64
    %xi = insertelement <$B x i64> poison, i64 %xz, i32 0
    %xs = shufflevector <$B x i64> %xi, <$B x i64> poison, <$B x i32> zeroinitializer
    %pr = mul <$B x i64> %mz, %xs
    %t = sub <$B x i64> %a, %pr
    %n = icmp slt <$B x i64> %t, zeroinitializer
    %qi = insertelement <$B x i64> poison, i64 %3, i32 0
    %qs = shufflevector <$B x i64> %qi, <$B x i64> poison, <$B x i32> zeroinitializer
    %ad = select <$B x i1> %n, <$B x i64> %qs, <$B x i64> zeroinitializer
    %r = add <$B x i64> %t, %ad
    store <$B x i64> %r, ptr %0, align 8
    ret void
    """
    quote
        Base.llvmcall(
            $ir,
            Cvoid,
            Tuple{Ptr{Int64}, NTuple{B, VecElement{UInt32}}, UInt32, Int64},
            pa,
            m,
            x,
            p2
        )
    end
end

# a[idx[k]][1..B] -= m[1..B] * c[k], kept in [0, p^2)
@inline function _blk_axpy!(
    A::Vector{Int64},
    T::Vector{UInt8},
    idx::Vector{I},
    c::Vector{C},
    m::NTuple{B, UInt32},
    p2::Int64
) where {I, C, B}
    mv = ntuple(l -> VecElement(m[l]), Val(B))
    n = length(idx)
    GC.@preserve A T idx c begin
        pa, pt, pi, pc = pointer(A), pointer(T), pointer(idx), pointer(c)
        for k in 1:n
            j = Int(unsafe_load(pi, k))
            unsafe_store!(pt, 0x01, j)
            _blk_vupd!(pa + 8 * (j - 1) * B, mv, unsafe_load(pc, k) % UInt32, p2)
        end
    end
    nothing
end

const _MI = Base.MultiplicativeInverses.UnsignedMultiplicativeInverse{UInt64}

@inline _blk_rem(x::Int64, mi::_MI) = (u = reinterpret(UInt64, x); u - div(u, mi) * mi.divisor)

@inline function _blk_mults(A::Vector{Int64}, base::Int, ::Val{B}, mi::_MI) where {B}
    @inbounds ntuple(l -> _blk_rem(A[base + l], mi) % UInt32, Val(B))
end

# first c' >= c with T[c'] != 0, or last + 1
@inline function _blk_next_touched(T::Vector{UInt8}, c::Int, last::Int)
    GC.@preserve T begin
        pt = pointer(T)
        @inbounds while c + 7 <= last
            w = unsafe_load(Ptr{UInt64}(pt + c - 1))
            w != 0 && return c + (trailing_zeros(w) >> 3)
            c += 8
        end
        @inbounds while c <= last
            T[c] != 0 && return c
            c += 1
        end
    end
    last + 1
end

# Lower part of the matrix reduced by blocks of B rows. learn: records the useful reducers and
# the rows not reduced to zero, as linalg_learn_reduce_matrix_lower_part!. apply: false as soon
# as a row reduces to zero, as linalg_apply_reduce_matrix_lower_part!.
# The dense block and its column marks, kept per trace from one call to the next: every call
# leaves them zero, but one that stops at a row reduced to zero (dirty)
mutable struct _BlkBuffers
    A::Vector{Int64}
    T::Vector{UInt8}
    dirty::Bool
end
const _BLK_BUFFERS = WeakKeyDict{Any, _BlkBuffers}()
const _BLK_LOCK = ReentrantLock()

function _blk_buffers(trace, ncols::Int, B::Int)
    buf = @lock _BLK_LOCK get!(() -> _BlkBuffers(Int64[], UInt8[], false), _BLK_BUFFERS, trace)
    if buf.dirty
        fill!(buf.A, 0)
        fill!(buf.T, 0x00)
        buf.dirty = false
    end
    if length(buf.A) < ncols * B
        n = length(buf.A)
        resize!(buf.A, ncols * B)
        fill!(view(buf.A, (n + 1):(ncols * B)), 0)
    end
    if length(buf.T) < ncols
        n = length(buf.T)
        resize!(buf.T, ncols)
        fill!(view(buf.T, (n + 1):ncols), 0x00)
    end
    buf
end

function linalg_reduce_matrix_lower_part_blocked!(
    trace::Union{Trace, Nothing},
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic,
    ::Val{B},
    learn::Bool
) where {CoeffType, B}
    _, ncols = size(matrix)
    nleft, _ = matrix_ncols_filled(matrix)
    _, nlow = matrix_nrows_filled(matrix)
    p = Int64(divisor(arithmetic))
    p2 = p * p
    mi = _MI(UInt64(p))

    pivots, row_idx_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # without a trace (the classic F4): rows reduced to zero allowed, nothing recorded, the
    # buffers kept with the basis
    rec = learn && trace !== nothing
    buf = _blk_buffers(trace === nothing ? basis : trace, ncols, B)
    A, T = buf.A, buf.T
    cand = Int[]
    vals = Int64[]
    not_reduced_to_zero = Int[]
    useful = falses(rec ? nleft : 0)
    umask = zeros(UInt32, rec ? nleft : 0)
    visited = Int[]

    @inbounds for r0 in 1:B:nlow
        nb = min(B, nlow - r0 + 1)
        lmin, from = ncols + 1, ncols + 1
        for l in 1:nb
            sup = matrix.lower_rows[r0 + l - 1]
            cf = basis.coeffs[row_idx_to_coeffs[r0 + l - 1]]
            for k in eachindex(sup)
                j = Int(sup[k])
                A[(j - 1) * B + l] = Int64(cf[k])
                T[j] = 0x01
            end
            lmin = min(lmin, Int(sup[1]))
            from = min(from, max(Int(sup[1]), nleft + 1))
        end
        # left block: every column has a reducer from the upper part
        c = lmin - 1
        while true
            c = _blk_next_touched(T, c + 1, nleft)
            c > nleft && break
            base = (c - 1) * B
            m = _blk_mults(A, base, Val(B), mi)
            mask = UInt32(0)
            for l in 1:B
                mask |= UInt32(!iszero(m[l])) << (l - 1)
            end
            if mask != 0
                _blk_axpy!(A, T, pivots[c], basis.coeffs[matrix.upper_to_coeffs[c]], m, p2)
                if rec
                    umask[c] = mask
                    push!(visited, c)
                end
            end
            for l in 1:B
                A[base + l] = 0
            end
            T[c] = 0x00
        end
        # pivots found by the previous blocks
        for c in from:ncols
            (T[c] == 0 || !isassigned(pivots, c)) && continue
            base = (c - 1) * B
            m = _blk_mults(A, base, Val(B), mi)
            if any(!iszero, m)
                _blk_axpy!(A, T, pivots[c], matrix.some_coeffs[matrix.lower_to_coeffs[c]], m, p2)
            end
            for l in 1:B
                A[base + l] = 0
            end
        end
        empty!(cand)
        for c in from:ncols
            T[c] != 0 && push!(cand, c)
        end
        # the rows of the block, one after the other
        nzmask = UInt32(0)
        for l in 1:nb
            sup = ColumnLabel[]
            empty!(vals)
            for c in cand
                v = Int64(_blk_rem(A[(c - 1) * B + l], mi))
                if !iszero(v)
                    push!(sup, ColumnLabel(c))
                    push!(vals, v)
                end
            end
            if isempty(sup)
                learn || (buf.dirty = true; return false)
                continue
            end
            nzmask |= UInt32(1) << (l - 1)
            inv = invmod(vals[1], p)
            cf = Vector{CoeffType}(undef, length(vals))
            for k in eachindex(vals)
                cf[k] = CoeffType(rem(vals[k] * inv, p))
            end
            i = r0 + l - 1
            learn && push!(not_reduced_to_zero, i)
            lead = Int(sup[1])
            matrix.some_coeffs[i] = cf
            pivots[lead] = sup
            matrix.lower_to_coeffs[lead] = i
            if l < nb
                base = (lead - 1) * B
                m = ntuple(q -> q <= l ? zero(UInt32) : _blk_rem(A[base + q], mi) % UInt32, Val(B))
                any(!iszero, m) && _blk_axpy!(A, T, sup, cf, m, p2)
            end
        end
        for c in cand
            base = (c - 1) * B
            for l in 1:B
                A[base + l] = 0
            end
            T[c] = 0x00
        end
        if rec
            for c in visited
                (umask[c] & nzmask) != 0 && (useful[c] = true)
                umask[c] = 0
            end
            empty!(visited)
        end
    end

    if rec
        useful_reducers_sorted = findall(useful)
        push!(
            trace.matrix_infos,
            (nup=matrix.nrows_filled_upper, nlow=matrix.nrows_filled_lower, ncols=ncols)
        )
        push!(trace.matrix_nonzeroed_rows, not_reduced_to_zero)
        push!(
            trace.matrix_upper_rows,
            (
                map(f -> matrix.upper_to_coeffs[f], useful_reducers_sorted),
                map(f -> matrix.upper_to_mult[f], useful_reducers_sorted)
            )
        )
        push!(
            trace.matrix_lower_rows,
            (row_idx_to_coeffs[not_reduced_to_zero], matrix.lower_to_mult[not_reduced_to_zero])
        )
    end
    true
end

# domain of the 32 x 32 -> 64 kernels: one prime p < 2^31, coefficients UInt32; the vector
# kernel is written with opaque pointers (LLVM 16, Julia 1.11 and later)
_fast_arith(a) =
    VERSION >= v"1.11" &&
    a isa AbstractArithmeticZp{UInt64, UInt32} &&
    divisor(a) < (UInt64(1) << 31)
# rows reduced together (2, 4, 8 or 16); 0 disables the blocked reduction
const LINALG_BLOCK = Ref(8)

_blk_ok(arithmetic) = LINALG_BLOCK[] > 0 && _fast_arith(arithmetic)

# below this many rows to reduce, the blocks cost more than they save
const BLOCK_MIN_ROWS = Ref(0)
_blk_use(arithmetic, matrix) = _blk_ok(arithmetic) && matrix.nrows_filled_lower >= BLOCK_MIN_ROWS[]

function _blk_dispatch(trace, matrix, basis, arithmetic, learn)
    b = LINALG_BLOCK[]
    b == 2 && return linalg_reduce_matrix_lower_part_blocked!(
        trace,
        matrix,
        basis,
        arithmetic,
        Val(2),
        learn
    )
    b == 4 && return linalg_reduce_matrix_lower_part_blocked!(
        trace,
        matrix,
        basis,
        arithmetic,
        Val(4),
        learn
    )
    b == 16 && return linalg_reduce_matrix_lower_part_blocked!(
        trace,
        matrix,
        basis,
        arithmetic,
        Val(16),
        learn
    )
    linalg_reduce_matrix_lower_part_blocked!(trace, matrix, basis, arithmetic, Val(8), learn)
end

# The randomized reduction of the lower part (linalg_randomized_reduce_matrix_lower_part!), B
# random combinations of the rows of a group at a time: loaded and reduced by the pivots
# together, then made pivots one after the other. A group is done when a combination
# reduces to zero, or when it has given as many pivots as it has rows.
function linalg_randomized_reduce_matrix_lower_part_blocked!(
    matrix::MacaulayMatrix{CoeffType},
    basis::Basis{CoeffType},
    arithmetic,
    rng::AbstractRNG,
    ::Val{B}
) where {CoeffType, B}
    _, ncols = size(matrix)
    nleft, _ = matrix_ncols_filled(matrix)
    _, nlow = matrix_nrows_filled(matrix)
    p = Int64(divisor(arithmetic))
    p2 = p * p
    mi = _MI(UInt64(p))

    pivots, row_idx_to_coeffs = linalg_prepare_matrix_pivots!(matrix)
    resize!(matrix.some_coeffs, nlow)

    # groups of RANDOM_GROUP rows, or as linalg_randomized_reduce_matrix_lower_part!
    if RANDOM_GROUP[] > 0
        rowsperblock = min(nlow, RANDOM_GROUP[])
        nblocks = cld(nlow, rowsperblock)
    else
        nblocks = linalg_nblocks_in_randomized(nlow)
        rowsperblock = div(nlow, nblocks) + (nlow % nblocks == 0 ? 0 : 1)
    end

    buf = _blk_buffers(basis, ncols, B)
    A, T = buf.A, buf.T
    cand = Int[]
    vals = Int64[]

    @inbounds for g in 1:nblocks
        first_row = (g - 1) * rowsperblock + 1
        last_row = min(nlow, g * rowsperblock)
        nrows = last_row - first_row + 1
        nrows <= 0 && continue
        newpivots, done = 0, false
        while !done && newpivots < nrows
            # B random combinations of the group rows: c_l += r * row, as c_l -= (p - r) * row
            lmin = ncols + 1
            for k in first_row:last_row
                sup = matrix.lower_rows[k]
                cf = basis.coeffs[row_idx_to_coeffs[k]]
                m = ntuple(l -> UInt32(p - 1 - rand(rng, 0:(p - 2))), Val(B))
                _blk_axpy!(A, T, sup, cf, m, p2)
                lmin = min(lmin, Int(sup[1]))
            end
            # left block: every column has a reducer from the upper part
            c = lmin - 1
            while true
                c = _blk_next_touched(T, c + 1, nleft)
                c > nleft && break
                base = (c - 1) * B
                m = _blk_mults(A, base, Val(B), mi)
                any(!iszero, m) &&
                    _blk_axpy!(A, T, pivots[c], basis.coeffs[matrix.upper_to_coeffs[c]], m, p2)
                for l in 1:B
                    A[base + l] = 0
                end
                T[c] = 0x00
            end
            # the pivots found so far in the lower part
            for c in (nleft + 1):ncols
                (T[c] == 0 || !isassigned(pivots, c)) && continue
                base = (c - 1) * B
                m = _blk_mults(A, base, Val(B), mi)
                any(!iszero, m) && _blk_axpy!(
                    A,
                    T,
                    pivots[c],
                    matrix.some_coeffs[matrix.lower_to_coeffs[c]],
                    m,
                    p2
                )
                for l in 1:B
                    A[base + l] = 0
                end
            end
            empty!(cand)
            for c in (nleft + 1):ncols
                T[c] != 0 && push!(cand, c)
            end
            # the combinations, one after the other
            for l in 1:B
                (done || newpivots >= nrows) && break
                sup = ColumnLabel[]
                empty!(vals)
                for c in cand
                    v = Int64(_blk_rem(A[(c - 1) * B + l], mi))
                    if !iszero(v)
                        push!(sup, ColumnLabel(c))
                        push!(vals, v)
                    end
                end
                if isempty(sup)
                    done = true
                    break
                end
                inv = invmod(vals[1], p)
                cf = Vector{CoeffType}(undef, length(vals))
                for k in eachindex(vals)
                    cf[k] = CoeffType(rem(vals[k] * inv, p))
                end
                newpivots += 1
                i = first_row - 1 + newpivots
                lead = Int(sup[1])
                matrix.some_coeffs[i] = cf
                pivots[lead] = sup
                matrix.lower_to_coeffs[lead] = i
                if l < B
                    base = (lead - 1) * B
                    m = ntuple(
                        q -> q <= l ? zero(UInt32) : _blk_rem(A[base + q], mi) % UInt32,
                        Val(B)
                    )
                    any(!iszero, m) && _blk_axpy!(A, T, sup, cf, m, p2)
                end
            end
            # the block left zero
            for c in cand
                base = (c - 1) * B
                for l in 1:B
                    A[base + l] = 0
                end
                T[c] = 0x00
            end
        end
    end
    true
end

# rows per group of the randomized reduction by blocks; 0: those of Groebner.jl
const RANDOM_GROUP = Ref(256)

function _blk_random_dispatch(matrix, basis, arithmetic, rng)
    b = LINALG_BLOCK[]
    b == 4 && return linalg_randomized_reduce_matrix_lower_part_blocked!(
        matrix,
        basis,
        arithmetic,
        rng,
        Val(4)
    )
    b == 16 && return linalg_randomized_reduce_matrix_lower_part_blocked!(
        matrix,
        basis,
        arithmetic,
        rng,
        Val(16)
    )
    linalg_randomized_reduce_matrix_lower_part_blocked!(matrix, basis, arithmetic, rng, Val(8))
end

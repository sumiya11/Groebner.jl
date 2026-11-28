# Replayed apply: the structure of every matrix is kept after a first apply, the next ones only
# redo the linear algebra and copy the coefficients.

# false disables the replay
const REPLAY = Ref(true)

mutable struct ReplayState
    steps::Vector{Any}
    auto::Any
    done::Bool
end

# one state per trace, held weakly; locked, traces may be applied in different threads
const _REPLAY = WeakKeyDict{Any, ReplayState}()
const _REPLAY_LOCK = ReentrantLock()

_replay_state(trace) = @lock _REPLAY_LOCK get(_REPLAY, trace, nothing)
_replay_forget!(trace) = (@lock _REPLAY_LOCK delete!(_REPLAY, trace); nothing)

function _replay_auto(trace)
    st = _replay_state(trace)
    st === nothing ? nothing : st.auto
end

struct ReplayStep
    upper_rows::Vector{Vector{ColumnLabel}}
    lower_rows::Vector{Vector{ColumnLabel}}
    upper_to_coeffs::Vector{Int}
    lower_to_coeffs::Vector{Int}
    ncols_left::Int
    ncols_right::Int
    out::Vector{Vector{ColumnLabel}}
end

# only during the recording apply
function _replay_capture_start(trace, matrix, arithmetic)
    (REPLAY[] && _fast_arith(arithmetic)) || return nothing
    st = _replay_state(trace)
    (st === nothing || st.done) && return nothing
    (
        map(copy, matrix.upper_rows[1:(matrix.nrows_filled_upper)]),
        map(copy, matrix.lower_rows[1:(matrix.nrows_filled_lower)]),
        copy(matrix.upper_to_coeffs[1:(matrix.nrows_filled_upper)]),
        copy(matrix.lower_to_coeffs[1:(matrix.nrows_filled_lower)]),
        matrix.ncols_left,
        matrix.ncols_right
    )
end

function _replay_capture_end!(trace, cap, matrix)
    out = map(copy, matrix.lower_rows[1:(matrix.npivots)])
    push!(_replay_state(trace).steps, ReplayStep(cap..., out))
    nothing
end

# the recorded steps, or nothing: then this apply records, the next ones replay
function _replay_get(trace, iters_total, arithmetic)
    (REPLAY[] && _fast_arith(arithmetic)) || return nothing
    st = _replay_state(trace)
    st !== nothing && st.done && return st.steps
    @lock _REPLAY_LOCK (_REPLAY[trace] = ReplayState(Any[], nothing, false))
    nothing
end

function _replay_mark_done!(trace)
    st = _replay_state(trace)
    st === nothing || (st.done = true)
    nothing
end

function _replay_iteration!(st::ReplayStep, trace, basis, matrix, params)
    nup, nlow = length(st.upper_rows), length(st.lower_rows)
    matrix.upper_rows = copy(st.upper_rows)
    matrix.lower_rows = copy(st.lower_rows)
    matrix.upper_to_coeffs = copy(st.upper_to_coeffs)
    matrix.lower_to_coeffs = copy(st.lower_to_coeffs)
    matrix.ncols_left = st.ncols_left
    matrix.ncols_right = st.ncols_right
    matrix.nrows_filled_upper = nup
    matrix.nrows_filled_lower = nlow
    flag = linalg_main!(matrix, basis, params, trace, linalg=LinearAlgebra(:apply, :sparse))
    !flag && return false
    matrix.npivots == length(st.out) || return false
    basis_resize_if_needed!(basis, matrix.npivots)
    crs = basis.n_processed
    @inbounds for i in 1:(matrix.npivots)
        r = matrix.lower_rows[i]
        r == st.out[i] || return false
        basis.coeffs[crs + i] = matrix.some_coeffs[matrix.lower_to_coeffs[r[1]]]
    end
    basis.n_filled += matrix.npivots
    true
end

# Final autoreduction replayed: only the output elements are reduced, by the unreduced
# reducers of the autoreduction matrix, column by column; the other rows are not touched.
struct ReplayAuto
    upper_rows::Vector{Vector{ColumnLabel}}
    upper_to_coeffs::Vector{Int}
    ncols::Int
    crs::Int
    npivots::Int
    outs::Vector{Tuple{Int, Vector{ColumnLabel}}}   # (index in the basis, support in columns)
end

function _replay_capture_auto!(trace, cap, matrix, crs)
    up, _, utc, _, nl, nr = cap
    outs = Tuple{Int, Vector{ColumnLabel}}[]
    for o in trace.output_nonredundant_indices
        push!(outs, (o, copy(matrix.lower_rows[o - crs])))
    end
    _replay_state(trace).auto = ReplayAuto(up, utc, nl + nr, crs, matrix.npivots, outs)
    nothing
end

function _replay_autoreduce!(ra::ReplayAuto, trace, basis::Basis{C}, params) where {C}
    B = 8
    ncols = ra.ncols
    p = Int64(divisor(params.arithmetic))
    p2 = p * p
    mi = _MI(UInt64(p))
    pivots = Vector{Vector{ColumnLabel}}(undef, ncols)
    colcoef = zeros(Int, ncols)
    @inbounds for k in eachindex(ra.upper_rows)
        r = ra.upper_rows[k]
        pivots[r[1]] = r
        colcoef[r[1]] = ra.upper_to_coeffs[k]
    end
    A = zeros(Int64, ncols * B)
    T = zeros(UInt8, ncols)
    outs = ra.outs
    nout = length(outs)
    sup = ColumnLabel[]
    vals = Int64[]
    @inbounds for r0 in 1:B:nout
        nb = min(B, nout - r0 + 1)
        cmin = ncols + 1
        leads = ntuple(l -> l <= nb ? Int(outs[r0 + l - 1][2][1]) : 0, Val(8))
        for l in 1:nb
            c = leads[l]
            cmin = min(cmin, c)
            r = pivots[c]
            cf = basis.coeffs[colcoef[c]]
            for k in 2:length(r)
                j = Int(r[k])
                A[(j - 1) * B + l] = Int64(cf[k])
                T[j] = 0x01
            end
        end
        c = cmin
        while true
            c = _blk_next_touched(T, c + 1, ncols)
            c > ncols && break
            isassigned(pivots, c) || continue
            base = (c - 1) * B
            m = _blk_mults(A, base, Val(B), mi)
            if any(!iszero, m)
                _blk_axpy!(A, T, pivots[c], basis.coeffs[colcoef[c]], m, p2)
            end
            for l in 1:B
                A[base + l] = 0
            end
            T[c] = 0x00
        end
        for l in 1:nb
            o, rec = outs[r0 + l - 1]
            lead = leads[l]
            invl = invmod(Int64(basis.coeffs[colcoef[lead]][1]), p)
            empty!(sup)
            empty!(vals)
            push!(sup, ColumnLabel(lead))
            push!(vals, 1)
            c = cmin
            while true
                c = _blk_next_touched(T, c + 1, ncols)
                c > ncols && break
                v = Int64(_blk_rem(A[(c - 1) * B + l], mi))
                if !iszero(v)
                    push!(sup, ColumnLabel(c))
                    push!(vals, rem(v * invl, p))
                end
            end
            sup == rec || return false
            basis.coeffs[o] = C.(vals)
        end
        c = cmin
        while true
            c = _blk_next_touched(T, c + 1, ncols)
            c > ncols && break
            base = (c - 1) * B
            for l in 1:B
                A[base + l] = 0
            end
            T[c] = 0x00
        end
    end
    basis.n_filled = ra.npivots + ra.crs
    basis.n_processed = ra.npivots
    true
end

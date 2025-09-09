# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Auxiliary functions wrapping libflint
# Adapted from Nemo.jl/src/flint/fmpq.jl

function flint_mpz_to_zz!(b::Nemo.ZZRingElem, a::BigInt)
    @ccall Nemo.libflint.fmpz_set_mpz(b::Ref{Nemo.ZZRingElem}, a::Ref{BigInt})::Nothing
    b
end

function flint_qq_to_mpq!(b::Rational{BigInt}, a::Nemo.QQFieldElem)
    @ccall Nemo.libflint.fmpq_get_mpz_frac(
        b.num::Ref{BigInt},
        b.den::Ref{BigInt},
        a::Ref{Nemo.QQFieldElem}
    )::Nothing
    b
end

function flint_reconstruct_fmpq_2!(
    c::Nemo.QQFieldElem,
    a::Nemo.ZZRingElem,
    m::Nemo.ZZRingElem,
    N::Nemo.ZZRingElem,
    D::Nemo.ZZRingElem
)
    success = Bool(
        @ccall Nemo.libflint.fmpq_reconstruct_fmpz_2(
            c::Ref{Nemo.QQFieldElem},
            a::Ref{Nemo.ZZRingElem},
            m::Ref{Nemo.ZZRingElem},
            N::Ref{Nemo.ZZRingElem},
            D::Ref{Nemo.ZZRingElem}
        )::Cint
    )
    success, c
end

###
# Rational reconstruction

# Returns the largest integer N (possibly off by 1) such that 2 N^2 < m.
function ratrec_reconstruction_bound(m::BigInt)
    isqrt((m >> UInt64(1)) - 1)
end

function ratrec_nemo!(
    res::Rational{BigInt},
    c::Nemo.QQFieldElem,
    a::Nemo.ZZRingElem,
    m::Nemo.ZZRingElem,
    N::Nemo.ZZRingElem,
    D::Nemo.ZZRingElem
)
    success, _ = flint_reconstruct_fmpq_2!(c, a, m, N, D)
    flint_qq_to_mpq!(res, c)
    success, res
end

function ratrec_nemo(a::Nemo.ZZRingElem, m::Nemo.ZZRingElem, N::Nemo.ZZRingElem, D::Nemo.ZZRingElem)
    ratrec_nemo!(zero(Rational{BigInt}), zero(Nemo.QQFieldElem), a, m, N, D)
end

###
# Element-wise rational reconstruction

# Table of rationals need not be initialized.
# Reconstructs only the witness set.
@timeit _TIMER function ratrec_vec_partial!(
    table_qq::Vector{Vector{Rational{BigInt}}},
    table_zz::Vector{Vector{BigInt}},
    modulo::BigInt,
    witness_set::Vector{Tuple{Int, Int}},
    mask::Vector{BitVector}
)
    @invariant length(table_qq) == length(table_zz)
    @invariant modulo > 1

    bound = ratrec_reconstruction_bound(modulo)
    nemo_bound = Nemo.ZZRingElem(bound)
    nemo_modulo = Nemo.ZZRingElem(modulo)
    rem_nemo = Nemo.ZZRingElem(0)
    nemo_buf = zero(Nemo.QQFieldElem)

    @inbounds for k in 1:length(witness_set)
        i, j = witness_set[k]
        @invariant 1 <= i <= length(table_zz) && 1 <= j <= length(table_zz[i])
        flint_mpz_to_zz!(rem_nemo, table_zz[i][j])
        success, _ =
            ratrec_nemo!(table_qq[i][j], nemo_buf, rem_nemo, nemo_modulo, nemo_bound, nemo_bound)
        !success && return false
        mask[i][j] = true
    end

    true
end

function _ratrec_vec_full!(
    table_qq::Vector{Vector{Rational{BigInt}}},
    table_zz::Vector{Vector{BigInt}},
    modulo::BigInt,
    mask::Vector{BitVector},
    rem_nemo::Nemo.ZZRingElem,
    nemo_buf::Nemo.QQFieldElem,
    nemo_modulo::Nemo.ZZRingElem,
    nemo_bound_N::Nemo.ZZRingElem,
    nemo_bound_D::Nemo.ZZRingElem,
    chunk::Vector{Int}
)
    @inbounds for i in chunk
        @invariant length(table_zz[i]) == length(table_qq[i])
        for j in 1:length(table_zz[i])
            mask[i][j] && continue
            flint_mpz_to_zz!(rem_nemo, table_zz[i][j])
            @invariant 0 <= rem_nemo < modulo
            success, _ = ratrec_nemo!(
                table_qq[i][j],
                nemo_buf,
                rem_nemo,
                nemo_modulo,
                nemo_bound_N,
                nemo_bound_D
            )
            !success && return false
        end
    end
    true
end

# Table of rationals need not be initialized.
@timeit _TIMER function ratrec_vec_full!(
    table_qq::Vector{Vector{Rational{BigInt}}},
    table_zz::Vector{Vector{BigInt}},
    modulo::BigInt,
    mask::Vector{BitVector};
    n_tasks=1
)
    @invariant length(table_qq) == length(table_zz)
    @invariant modulo > 1

    bound = ratrec_reconstruction_bound(modulo)
    nemo_bound = Nemo.ZZRingElem(bound)
    nemo_modulo = Nemo.ZZRingElem(modulo)
    # rem_nemo = Nemo.ZZRingElem(0)

    n_tasks = min(n_tasks, length(table_zz))
    chunk_size = max(1, div(length(table_zz), n_tasks))
    data_chunks = [
        [i + n_tasks * (j - 1) for j in 1:chunk_size if i + n_tasks * (j - 1) <= length(table_zz)] for i in 1:n_tasks
    ]

    tasks = Vector{Task}(undef, length(data_chunks))
    for (tid, chunk) in enumerate(data_chunks)
        task = @spawn begin
            local rem_nemo = zero(Nemo.ZZRingElem)
            local nemo_buf = zero(Nemo.QQFieldElem)
            _ratrec_vec_full!(
                table_qq,
                table_zz,
                modulo,
                mask,
                rem_nemo,
                nemo_buf,
                nemo_modulo,
                nemo_bound,
                nemo_bound,
                chunk
            )
        end
        tasks[tid] = task
    end
    all(map(fetch, tasks))
    # for task in tasks
    #     wait(task)
    # end
    # nothing

    # @inbounds for i in 1:length(table_zz)
    #     @invariant length(table_zz[i]) == length(table_qq[i])
    #     for j in 1:length(table_zz[i])
    #         mask[i][j] && continue
    #         mpz_to_zz!(rem_nemo, table_zz[i][j])
    #         @invariant 0 <= rem_nemo < modulo
    #         success, _ = ratrec_nemo!(table_qq[i][j], rem_nemo, nemo_modulo, nemo_bound, nemo_bound)
    #         !success && return false
    #     end
    # end

    # true
end

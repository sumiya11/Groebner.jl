# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Linear algebra, main entry points

"""
    linalg_main!

Given a `MacaulayMatrix` of the form

    | A B |
    | C D |

linearly reduces the `CD` part with respect to the known pivots from the `AB`
part, and interreduces (or, autoreduces) the result inplace.

In other words, computes the reduced row echelon form of

    D - C inv(A) B.

Returns `true` if successful and `false` otherwise.
"""
function linalg_main!(
    matrix::MacaulayMatrix,
    basis::Basis,
    params::AlgorithmParameters,
    trace=nothing;
    linalg=nothing
)
    @invariant matrix_well_formed(matrix)

    rng = params.rng
    arithmetic = params.arithmetic
    if isnothing(linalg)
        linalg = params.linalg
    end
    changematrix = params.changematrix

    # Decide if multi-threading should be used. Further in dispatch, the backend
    # may opt to NOT use multi-threading regardless of this setting.
    threaded = if params.threaded_f4 === :yes && nthreads() > 1
        # use multi-threading if explicitly requested
        :yes
    elseif params.threaded_multimodular === :yes && linalg.algorithm === :learn && nthreads() > 1
        # TODO: we do not use threading at the learn stage, but we could. One of
        # the obstactles here are frequent allocations and GC.
        :no
    elseif params.threaded_f4 === :auto && nthreads() > 1
        # try to use multi-threading by default whenever possible
        :yes
    else
        :no
    end

    flag = if !isnothing(trace)
        _linalg_main_with_trace!(trace, matrix, basis, linalg, threaded, changematrix, arithmetic, rng)
    else
        _linalg_main!(matrix, basis, linalg, threaded, changematrix, arithmetic, rng)
    end

    flag
end

"""
    linalg_autoreduce!

Interreduces the rows of the given `MacaulayMatrix`.

Returns `true` if successful and `false` otherwise.
"""
function linalg_autoreduce!(
    matrix::MacaulayMatrix,
    basis::Basis,
    params::AlgorithmParameters,
    trace=nothing;
    linalg=nothing
)
    @invariant matrix_well_formed(matrix)

    arithmetic = params.arithmetic
    if isnothing(linalg)
        linalg = params.linalg
    end
    changematrix = params.changematrix

    flag = if !isnothing(trace)
        _linalg_autoreduce_with_trace!(trace, matrix, basis, linalg, changematrix, arithmetic)
    else
        _linalg_autoreduce!(matrix, basis, linalg, changematrix, arithmetic)
    end

    flag
end

"""
    linalg_normalform!

Given a `MacaulayMatrix` of the form

    | A B |
    | C D |

linearly reduces the `CD` part with respect to the pivots in the `AB` part.

In contrast to `linalg_main!`, this function does not perform the final
interreduction (or, autoreduction) of the matrix rows.

Returns `true` if successful and `false` otherwise.
"""
function linalg_normalform!(matrix::MacaulayMatrix, basis::Basis, arithmetic::AbstractArithmetic)
    @invariant matrix_well_formed(matrix)
    _linalg_normalform!(matrix, basis, arithmetic)
end

"""
    linalg_isgroebner!

Given a `MacaulayMatrix` of the form

    | A B |
    | C D |

returns `true` if

    D - C inv(A) B

is zero and `false`, otherwise.
"""
function linalg_isgroebner!(matrix::MacaulayMatrix, basis::Basis, arithmetic::AbstractArithmetic)
    @invariant matrix_well_formed(matrix)
    _linalg_isgroebner!(matrix, basis, arithmetic)
end

###
# Further dispatch between linear algebra backends

function linalg_should_use_threading(matrix, linalg, threaded)
    @assert threaded !== :auto
    nup, nlow, nleft, nright = matrix_block_sizes(matrix)
    if threaded === :yes
        # Opt for single thread for small matrices
        if nlow < 2 * nthreads()
            false
        elseif nleft + nright < 1_000
            false
        else
            true
        end
    else
        false
    end
end

function _linalg_main!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    threaded::Symbol,
    changematrix::Bool,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    flag = if changematrix
        linalg_deterministic_sparse_changematrix!(matrix, basis, linalg, arithmetic)
    elseif linalg.algorithm === :deterministic
        if linalg_should_use_threading(matrix, linalg, threaded)
            linalg_deterministic_sparse_threaded!(matrix, basis, linalg, arithmetic)
        else
            linalg_deterministic_sparse!(matrix, basis, linalg, arithmetic)
        end
    elseif linalg.algorithm === :randomized
        if linalg_should_use_threading(matrix, linalg, threaded)
            linalg_randomized_sparse_threaded!(matrix, basis, linalg, arithmetic, rng)
        else
            linalg_randomized_sparse!(matrix, basis, linalg, arithmetic, rng)
        end
    else
        throw(DomainError("Cannot pick a suitable linear algebra backend"))
    end

    flag
end

function _linalg_main_with_trace!(
    trace::Trace,
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    threaded::Symbol,
    changematrix::Bool,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    flag = if linalg.algorithm === :learn
        # At the moment, this never opts for multi-threading
        if false && threaded === :yes
            linalg_learn_sparse_threaded!(trace, matrix, basis, arithmetic)
        else
            linalg_learn_sparse!(trace, matrix, basis, arithmetic)
        end
    else
        @assert linalg.algorithm === :apply
        linalg_apply_sparse!(trace, matrix, basis, arithmetic)
    end

    flag
end

function _linalg_autoreduce!(
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    changematrix::Bool,
    arithmetic::AbstractArithmetic
)
    if changematrix
        linalg_deterministic_sparse_interreduction_changematrix!(matrix, basis, arithmetic)
    else
        linalg_deterministic_sparse_interreduction!(matrix, basis, arithmetic)
    end
end

function _linalg_autoreduce_with_trace!(
    trace::Trace,
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg::LinearAlgebra,
    changematrix::Bool,
    arithmetic::AbstractArithmetic
)
    @assert !changematrix
    sort_matrix_upper_rows!(matrix)

    if linalg.algorithm === :learn
        linalg_learn_deterministic_sparse_interreduction!(trace, matrix, basis, arithmetic)
    else
        @assert linalg.algorithm === :apply
        linalg_apply_deterministic_sparse_interreduction!(trace, matrix, basis, arithmetic)
    end

    true
end

function _linalg_normalform!(matrix::MacaulayMatrix, basis::Basis, arithmetic::AbstractArithmetic)
    sort_matrix_upper_rows!(matrix)
    linalg_reduce_matrix_lower_part_invariant_pivots!(matrix, basis, arithmetic)
end

function _linalg_isgroebner!(matrix::MacaulayMatrix, basis::Basis, arithmetic::AbstractArithmetic)
    sort_matrix_upper_rows!(matrix)
    sort_matrix_lower_rows!(matrix)
    linalg_reduce_matrix_lower_part_all_zero!(matrix, basis, arithmetic)
end

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
    tasks = params.tasks
    flag = if !isnothing(trace)
        _linalg_main_with_trace!(trace, matrix, basis, linalg, tasks, changematrix, arithmetic, rng)
    else
        _linalg_main!(matrix, basis, linalg, tasks, changematrix, arithmetic, rng)
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
function linalg_isgroebner!(matrix::MacaulayMatrix, basis::Basis, params::AlgorithmParameters)
    @invariant matrix_well_formed(matrix)
    flag = if params.linalg.algorithm === :deterministic
        _linalg_isgroebner_deterministic!(matrix, basis, params.arithmetic)
    elseif params.linalg.algorithm === :randomized
        _linalg_isgroebner_randomized!(matrix, basis, params.arithmetic, params.rng)
    else
        throw(DomainError("Cannot pick a suitable linear algebra backend"))
    end
    flag
end

###
# Further dispatch between linear algebra backends

function linalg_should_use_threading(matrix, linalg, tasks)
    nup, nlow, nleft, nright = matrix_block_sizes(matrix)
    if tasks > 1
        # Opt for single thread for small matrices
        if nlow < 2 * tasks
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
    tasks::Int,
    changematrix::Bool,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    flag = if changematrix
        linalg_deterministic_sparse_changematrix!(matrix, basis, linalg, arithmetic)
    elseif linalg.algorithm === :deterministic
        if linalg_should_use_threading(matrix, linalg, tasks)
            linalg_deterministic_sparse_threaded!(matrix, basis, linalg, tasks, arithmetic)
        else
            linalg_deterministic_sparse!(matrix, basis, linalg, arithmetic)
        end
    elseif linalg.algorithm === :randomized
        if linalg_should_use_threading(matrix, linalg, tasks)
            linalg_randomized_sparse_threaded!(matrix, basis, linalg, tasks, arithmetic, rng)
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
    tasks::Int,
    changematrix::Bool,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    flag = if linalg.algorithm === :learn
        # At the moment, this never opts for multi-threading
        if linalg_should_use_threading(matrix, linalg, tasks)
            linalg_learn_sparse_threaded!(trace, matrix, basis, tasks, arithmetic)
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
    linalg_reduce_matrix_lower_part_do_not_modify_pivots!(matrix, basis, arithmetic)
end

function _linalg_isgroebner_deterministic!(
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic
)
    sort_matrix_upper_rows!(matrix)
    sort_matrix_lower_rows!(matrix)
    linalg_reduce_matrix_lower_part_all_zero!(matrix, basis, arithmetic)
end

function _linalg_isgroebner_randomized!(
    matrix::MacaulayMatrix,
    basis::Basis,
    arithmetic::AbstractArithmetic,
    rng::AbstractRNG
)
    sort_matrix_upper_rows!(matrix)
    sort_matrix_lower_rows!(matrix)
    linalg_randomized_reduce_matrix_lower_part_all_zero!(matrix, basis, arithmetic, rng)
end

# Macaulay matrix

# MacaulayMatrix is a sparse matrix with columns labeled by distinct monomials.
# MacaulayMatrix represents a block matrix of the following structure:
#
#   | A  B |
#   | C  D | 
#
# The sparsity structure of the matrix looks similar to:
#
#           A          B
#     ............ | ......
#       .......... | ......
#         ........ | ......
#  A        ...... | ......
#             .... |  .....
#               .. |   ....
#     ---------------------
#            ..... | ......
#  C      ........ | ......
#      ........... | ......
#
# NOTE: check out `repr_matrix` for printing a nice matrix representation.

# The primary action on a MacaulayMatrix is getting it into row echelon form.
#
# Usually, by construction,
# - Block A is already in row-echelon-form.
# - Permuting columns is not allowed.
#
# Upon completion, the matrix may look similar to:
#
#           A          B
#     ............ | ......
#       .......... | ......
#         ........ | ......
#  A        ...... | ......
#             .... |  .....
#               .. |   ....
#     ---------------------
#                  | ......
#  C               |   ....
#                  |     ..
#                  |      .
#

const ColumnLabel = Int32

# Tags for columns from different blocks
const NON_PIVOT_COLUMN = 0
const UNKNOWN_PIVOT_COLUMN = 1
const PIVOT_COLUMN = 2

mutable struct MacaulayMatrix{T <: Coeff}
    # A row in the matrix is represented sparsely with a vector of indices of
    # nonzero columns.
    # The coefficients of a matrix row can be stored in two different places:
    # - Either as a pointer to a polynomial in the basis that has the same
    #   coefficients.
    # - Or stored explicitly in the matrix.

    # Rows from the block AB.
    upper_rows::Vector{Vector{ColumnLabel}}
    # Rows from the block CD
    lower_rows::Vector{Vector{ColumnLabel}}

    # Explicitly stored coefficients of *some* rows
    coeffs::Vector{Vector{T}}

    # Maps column index in the range {1 ... ncols} to a corresponding monomial
    column_to_monom::Vector{MonomIdx}

    # The total number of allocated rows
    size::Int
    # The number of filled rows
    nrows::Int
    # The number of columns
    ncolumns::Int
    # The number of pivots, i.e. the number of linearly independent elements in
    # the block CD after the reduction by the block AB
    npivots::Int

    # The number of upper rows (in the AB block)
    nupper::Int
    # The number of lower rows (in the CD block)
    nlower::Int
    # The number of left columns (in the A and C blocks)
    nleft::Int
    # The number of right columns (in the B and D blocks)
    nright::Int

    # Maps column index to a pivot row.
    # pivots[i] may be undefined. In that case, there is no pivot row with the
    # leading column i
    pivots::Vector{Vector{ColumnLabel}}

    # Index of the row --> position of the array of coefficients for this row
    upper_to_coeffs::Vector{Int}
    lower_to_coeffs::Vector{Int}

    # Index of the row --> monomial multiplier
    upper_to_mult::Vector{MonomIdx}
    lower_to_mult::Vector{MonomIdx}
end

###
# MacaulayMatrix utilities

# Initializes an empty matrix with coefficients of type T
@timeit function initialize_matrix(ring::PolyRing, ::Type{T}) where {T <: Coeff}
    upper_rows = Vector{Vector{ColumnLabel}}(undef, 0)
    lower_rows = Vector{Vector{ColumnLabel}}(undef, 0)
    column_to_monom = Vector{MonomIdx}(undef, 0)
    coeffs = Vector{Vector{T}}(undef, 0)
    size = 0
    npivots = 0
    nrows = 0
    ncolumns = 0
    nupper = 0
    nlower = 0
    nleft = 0
    nright = 0
    pivots = Vector{Vector{ColumnLabel}}(undef, 0)
    upper_to_coeffs = Vector{Int}(undef, 0)
    lower_to_coeffs = Vector{Int}(undef, 0)
    upper_to_mult = Vector{MonomIdx}(undef, 0)
    lower_to_mult = Vector{MonomIdx}(undef, 0)
    MacaulayMatrix(
        upper_rows,
        lower_rows,
        coeffs,
        column_to_monom,
        size,
        nrows,
        ncolumns,
        npivots,
        nupper,
        nlower,
        nleft,
        nright,
        pivots,
        upper_to_coeffs,
        lower_to_coeffs,
        upper_to_mult,
        lower_to_mult
    )
end

# Returns a string with human-readable information about the matrix. Currently,
# used in logging (call, e.g., `groebner(system, loglevel=-3)` to check out)
function repr_matrix(matrix::MacaulayMatrix{T}) where {T}
    m, n = matrix.nrows, matrix.ncolumns
    m_A, n_A = matrix.nupper, matrix.nleft
    m_B, n_B = matrix.nupper, matrix.nright
    m_C, n_C = matrix.nlower, matrix.nleft
    m_D, n_D = matrix.nlower, matrix.nright
    nnz_A, nnz_B, nnz_C, nnz_D = 0, 0, 0, 0
    A_ref, A_rref = true, true
    # NOTE: probably want to adjust this when the terminal shrinks
    max_canvas_width = 40
    canvas = CanvasMatrix2x2(m_A, m_C, n_A, n_B, max_width=max_canvas_width)
    for i in 1:(matrix.nupper)
        row = matrix.upper_rows[i]
        if length(row) > 1
            A_rref = false
        end
        if row[1] < i
            A_ref = true
        end
        for j in 1:length(row)
            point!(canvas, i, row[j])
            if row[j] <= matrix.nleft
                nnz_A += 1
            else
                nnz_B += 1
            end
        end
    end
    for i in 1:(matrix.nlower)
        row = matrix.lower_rows[i]
        for j in 1:length(row)
            point!(canvas, matrix.nupper + i, row[j])
            if row[j] <= matrix.nleft
                nnz_C += 1
            else
                nnz_D += 1
            end
        end
    end
    nnz = nnz_A + nnz_B + nnz_C + nnz_D
    percent(x) = round(100 * x, digits=2)
    s = """
    $(typeof(matrix))
    $m x $n with $nnz nonzeros ($(percent(nnz / (m * n))) %)
    A: $(m_A) x $(n_A) with $(nnz_A) nonzeros (REF: $(A_ref), RREF: $(A_rref))
    B: $(m_B) x $(n_B) with $(nnz_B) nonzeros
    C: $(m_C) x $(n_C) with $(nnz_C) nonzeros
    D: $(m_D) x $(n_D) with $(nnz_D) nonzeros ($(percent(nnz_D / (m_D * n_D))) %)

    Sparsity pattern:

    $(canvas)
    """
    s
end

# Checks that the matrix is well formed.
# Generally this should be called on the entry to a linear algebra routine.
function matrix_well_formed(key, matrix::MacaulayMatrix)
    if !(matrix.size == length(matrix.upper_rows))
        return false
    end
    if !(matrix.ncolumns == matrix.nleft + matrix.nright)
        return false
    end
    if !(matrix.nrows == matrix.nupperper + matrix.nlower)
        return false
    end
    if key in (:in_reduction_apply!, :in_reduction!, :in_reduction_learn!)
        for i in 1:length(matrix.upper_rows)
            !isassigned(matrix.upper_rows, i) && continue
            if isempty(matrix.upper_rows[i])
                return false
            end
        end
        for i in 1:length(matrix.lower_rows)
            !isassigned(matrix.lower_rows, i) && continue
            if isempty(matrix.lower_rows[i])
                return false
            end
        end
    end
    true
end

# Refresh and initialize matrix
function reinitialize_matrix!(matrix::MacaulayMatrix{T}, size::Int) where {T}
    # TODO: do not reinitialize when the present size is enough 
    resize!(matrix.upper_rows, size * 2)
    resize!(matrix.lower_rows, size * 2)
    resize!(matrix.upper_to_coeffs, size * 2)
    resize!(matrix.lower_to_coeffs, size * 2)
    resize!(matrix.upper_to_mult, size * 2)
    resize!(matrix.lower_to_mult, size * 2)
    matrix.size = 2 * size
    matrix.ncolumns = 0
    matrix.nleft = 0
    matrix.nright = 0
    matrix.nupper = 0
    matrix.nlower = 0
    matrix
end

function resize_matrix_upper_part_if_needed!(matrix::MacaulayMatrix, to_add::Int)
    if matrix.size <= to_add
        while matrix.size <= to_add
            matrix.size *= 2
        end
        resize!(matrix.upper_rows, matrix.size)
        resize!(matrix.upper_to_coeffs, matrix.size)
        resize!(matrix.upper_to_mult, matrix.size)
    end
    nothing
end

###
# MacaulayMatrix linear algebra, high level

"""
    linear_algebra!(matrix, basis, params, [graph])

Produces a REF of the `matrix`. Returns `true` if successful.

Used in the F4 reduction.
"""
function linear_algebra! end

"""
    linear_algebra_reducegb!(matrix, basis, params, [graph])

Produces the RREF of the `matrix`. Returns `true` if successful.

Used in the autoreduction of the basis.
"""
function linear_algebra_reducegb! end

"""
    linear_algebra_normalform!(matrix, basis, params, [graph])

Puts the `matrix` in RREF, but does **not use** the rows from lower part of the
matrix as pivots for in reduction. Returns `true` if successful.

Used in the computation of normal forms.
"""
function linear_algebra_normalform! end

"""
    linear_algebra_isgroebner!(matrix, basis, params, [graph])

Checks that the lower part of the `matrix` reduces to zero w.r.t. the upper
part.

Used for assessing if all S-polynomials reduce to zero modulo a basis.
"""
function linear_algebra_isgroebner! end

function linear_algebra!(matrix, basis, params, graph=nothing; linalg::Symbol=:auto)
    # @invariant matrix_well_formed(:linear_algebra!, matrix)
    if linalg === :auto
        linalg = params.linalg
    end
    rng = params.rng
    arithmetic = params.arithmetic
    flag = if !isnothing(graph)
        linear_algebra!(graph, matrix, basis, linalg, arithmetic, rng)
    else
        linear_algebra!(matrix, basis, linalg, arithmetic, rng)
    end
    return flag
end

function linear_algebra!(matrix, basis, linalg, arithmetic, rng)
    resize!(matrix.coeffs, matrix.nlower)
    # @invariant matrix_well_formed(:none, matrix)
    flag = if linalg === :deterministic
        deterministic_sparse_linear_algebra!(matrix, basis, arithmetic)
    else
        @assert linalg === :randomized
        randomized_sparse_linear_algebra!(matrix, basis, arithmetic, rng)
    end
    flag
end

function linear_algebra!(
    graph::ComputationGraphF4,
    matrix::MacaulayMatrix,
    basis::Basis,
    linalg,
    arithmetic,
    rng
)
    resize!(matrix.coeffs, matrix.nlower)
    # @invariant matrix_well_formed(:none, matrix)
    flag = if linalg === :learn
        learn_sparse_linear_algebra!(graph, matrix, basis, arithmetic)
    else
        @assert linalg === :apply
        apply_sparse_linear_algebra!(graph, matrix, basis, arithmetic)
    end
    flag
end

function linear_algebra_reducegb!(
    matrix,
    basis,
    params,
    graph=nothing;
    linalg::Symbol=:auto
)
    sort_matrix_upper_rows!(matrix)

    arithmetic = params.arithmetic
    if linalg === :auto
        linalg = params.linalg
    end

    flag = if !isnothing(graph)
        linear_algebra_reducegb!(graph, matrix, basis, linalg, arithmetic)
    else
        linear_algebra_reducegb!(matrix, basis, linalg, arithmetic)
    end
    flag
end

function linear_algebra_reducegb!(
    matrix,
    basis,
    linalg,
    arithmetic::A
) where {A <: AbstractArithmetic}
    deterministic_sparse_interreduction!(matrix, basis, arithmetic)
end

function linear_algebra_reducegb!(
    graph,
    matrix,
    basis,
    linalg,
    arithmetic::A
) where {A <: AbstractArithmetic}
    if linalg === :learn
        learn_deterministic_sparse_interreduction!(graph, matrix, basis, arithmetic)
    else
        apply_deterministic_sparse_interreduction!(graph, matrix, basis, arithmetic)
    end
end

function linear_algebra_normalform!(matrix, basis, arithmetic)
    sort_matrix_upper_rows!(matrix)
    @log level = -3 "linear_algebra_normalform!"
    @log level = -3 repr_matrix(matrix)
    resize!(matrix.coeffs, matrix.nlower)
    reduce_matrix_lower_part_invariant_pivots!(matrix, basis, arithmetic)
end

function linear_algebra_isgroebner!(matrix, basis, arithmetic)
    sort_matrix_upper_rows!(matrix)
    sort_matrix_lower_rows!(matrix)
    @log level = -3 "linear_algebra_isgroebner!"
    @log level = -3 repr_matrix(matrix)
    reduce_matrix_lower_part_any_nonzero!(matrix, basis, arithmetic)
end

###
# MacaulayMatrix linear algebra, middle level

function deterministic_sparse_linear_algebra!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "deterministic_sparse_linear_algebra!"
    @log level = -3 repr_matrix(matrix)
    # Reduce CD with AB
    reduce_matrix_lower_part!(matrix, basis, arithmetic)
    # Interreduce CD
    interreduce_matrix_pivots!(matrix, basis, arithmetic)
end

function randomized_sparse_linear_algebra!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A,
    rng
) where {C <: Coeff, A <: AbstractArithmetic}
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "randomized_sparse_linear_algebra!"
    @log level = -3 repr_matrix(matrix)
    # Reduce CD with AB
    randomized_reduce_matrix_lower_part!(matrix, basis, arithmetic, rng)
    # Interreduce CD
    interreduce_matrix_pivots!(matrix, basis, arithmetic)
end

function learn_sparse_linear_algebra!(
    graph,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    sort_matrix_upper_rows!(matrix) # for the AB part
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "learn_sparse_linear_algebra!"
    @log level = -3 repr_matrix(matrix)
    # Reduce CD with AB
    learn_reduce_matrix_lower_part!(graph, matrix, basis, arithmetic)
    # Interreduce CD
    interreduce_matrix_pivots!(matrix, basis, arithmetic)
end

function apply_sparse_linear_algebra!(
    graph,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    # NOTE: here, we do not need to sort the rows, as they have already been
    # collected in the right order
    sort_matrix_lower_rows!(matrix) # for the CD part
    @log level = -3 "apply_sparse_linear_algebra!"
    @log level = -3 repr_matrix(matrix)
    # Reduce CD with AB
    flag = apply_reduce_matrix_lower_part!(graph, matrix, basis, arithmetic)
    if !flag
        return flag
    end
    # Interreduce CD
    apply_interreduce_matrix_pivots!(graph, matrix, basis, arithmetic)
end

function deterministic_sparse_interreduction!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    @log level = -3 "deterministic_sparse_interreduction!"
    @log level = -3 repr_matrix(matrix)
    # Prepare the matrix
    absolute_index_pivots_in_interreduction!(matrix, basis)
    # Interreduce AB
    interreduce_matrix_pivots!(matrix, basis, arithmetic, reversed_rows=true)
    true
end

function learn_deterministic_sparse_interreduction!(
    graph,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    @log level = -3 "learn_deterministic_sparse_interreduction!"
    @log level = -3 repr_matrix(matrix)
    # Prepare the matrix
    absolute_index_pivots_in_interreduction!(matrix, basis)
    # Interreduce AB
    learn_interreduce_matrix_pivots!(graph, matrix, basis, arithmetic, reversed_rows=true)
    true
end

function apply_deterministic_sparse_interreduction!(
    graph,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    @log level = -3 "apply_deterministic_sparse_interreduction!"
    @log level = -3 repr_matrix(matrix)
    # Prepare the matrix
    absolute_index_pivots_in_interreduction!(matrix, basis)
    # Interreduce AB
    flag = apply_interreduce_matrix_pivots!(
        graph,
        matrix,
        basis,
        arithmetic,
        reversed_rows=true
    )
    flag
end

###
# MacaulayMatrix linear algebra, low level

# Given a matrix of form 
#   A B
#   C D
# reduces the lower part CD with respect to the upper part AB.
# As a result, the matrix of the following form is produced:
#   A B
#   0 D' 
@timeit function reduce_matrix_lower_part!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    # Prepare the matrix
    pivots, row_index_to_coeffs = absolute_index_pivots!(matrix)
    ncols = matrix.ncolumns
    nlow = matrix.nlower
    pivots = matrix.pivots
    row = zeros(C, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(C)
    @inbounds for i in 1:nlow
        # Select the row from the lower part of the matrix to be reduced
        nnz_column_indices = matrix.lower_rows[i]
        # Locate the array of coefficients of this row.
        # NOTE: no copy of coefficients is needed
        nnz_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        # Load coefficients into a dense array
        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)
        # Reduce the row with respect to the known `pivots` from the upper part
        # of the matrix.
        # NOTE: this also does partial interreduction of the lower matrix rows.
        # Upon discovery of a new pivot from the lower part of the matrix, we
        # add it the pivot to the pool available of pivots
        first_nnz_column = nnz_column_indices[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            arithmetic,
            tmp_pos=-1
        )
        # If the row is fully reduced
        zeroed && continue

        # Store the new row in the matrix,
        # AND add it to the active pivots
        matrix.coeffs[i] = new_coeffs
        pivots[new_column_indices[1]] = new_column_indices
        # Set a reference to the coefficients of this row in the matrix
        matrix.lower_to_coeffs[new_column_indices[1]] = i

        normalize_sparse_row!(matrix.coeffs[i], arithmetic)

        new_column_indices, new_coeffs = new_empty_sparse_row(C)
    end
    true
end

# Given a matrix of form 
#   A B
#   C D
# reduces the lower part CD with respect to the upper part AB.
# As a result, the matrix of the following form is produced:
#   A B
#   0 D' 
#
# NOTE: Returns `false` if any row reduced to zero
function apply_reduce_matrix_lower_part!(
    graph,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    # Prepare the matrix
    pivots, row_index_to_coeffs = absolute_index_pivots!(matrix)
    ncols = matrix.ncolumns
    nlow = matrix.nlower
    pivots = matrix.pivots
    row = zeros(C, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(C)
    @inbounds for i in 1:nlow
        # Select the row from the lower part of the matrix to be reduced
        nnz_column_indices = matrix.lower_rows[i]
        # Locate the array of coefficients of this row.
        # NOTE: no copy of coefficients is needed
        nnz_coeffs = basis.coeffs[row_index_to_coeffs[i]]
        # Load coefficients into a dense array
        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)
        # Reduce the row with respect to the known `pivots` from the upper part
        # of the matrix.
        # NOTE: this also does partial interreduction of the lower matrix rows.
        # Upon discovery of a new pivot from the lower part of the matrix, we
        # add it the pivot to the pool available of pivots
        first_nnz_column = nnz_column_indices[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            arithmetic,
            tmp_pos=-1
        )
        # If the row is fully reduced
        if zeroed
            return false
        end

        # Store the new row in the matrix,
        # AND add it to the active pivots
        matrix.coeffs[i] = new_coeffs
        pivots[new_column_indices[1]] = new_column_indices
        # Set a reference to the coefficients of this row in the matrix
        matrix.lower_to_coeffs[new_column_indices[1]] = i

        normalize_sparse_row!(matrix.coeffs[i], arithmetic)

        new_column_indices, new_coeffs = new_empty_sparse_row(C)
    end
    true
end

# Interreduces the matrix.pivots 
@timeit function interreduce_matrix_pivots!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A;
    reversed_rows::Bool=false
) where {C <: Coeff, A <: AbstractArithmetic}
    pivots = matrix.pivots
    new_pivots = 0
    # Indices of rows that did no reduce to zero
    not_reduced_to_zero = Vector{Int}(undef, matrix.nright)
    row = zeros(C, matrix.ncolumns)
    resize!(matrix.lower_rows, matrix.nright)
    any_zeroed = false
    # for each column in the block D..
    @inbounds for i in 1:(matrix.nright)
        abs_column_idx = matrix.ncolumns - i + 1
        # Check if there is a row that starts at `abs_column_idx`
        !isassigned(pivots, abs_column_idx) && continue

        nnz_column_indices = pivots[abs_column_idx]
        if abs_column_idx <= matrix.nleft
            # upper part of matrix
            nnz_coeffs = basis.coeffs[matrix.upper_to_coeffs[abs_column_idx]]
        else
            # lower part of matrix
            nnz_coeffs = matrix.coeffs[matrix.lower_to_coeffs[abs_column_idx]]
        end
        @assert length(nnz_column_indices) == length(nnz_coeffs)
        load_sparse_row!(row, nnz_column_indices, nnz_coeffs)

        new_column_indices, new_coeffs = new_empty_sparse_row(C)
        first_nnz_column = nnz_column_indices[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            arithmetic,
            tmp_pos=first_nnz_column
        )
        # If the row reduced to zero
        if zeroed
            any_zeroed = true
            continue
        end

        new_pivots += 1
        not_reduced_to_zero[new_pivots] = i
        # update row and coeffs
        if !reversed_rows
            matrix.lower_rows[new_pivots] = new_column_indices
            matrix.coeffs[matrix.lower_to_coeffs[abs_column_idx]] = new_coeffs
            pivots[abs_column_idx] = matrix.lower_rows[new_pivots]
        else
            matrix.lower_rows[matrix.nrows - new_pivots + 1] = new_column_indices
            matrix.coeffs[matrix.lower_to_coeffs[abs_column_idx]] = new_coeffs
            pivots[abs_column_idx] = matrix.lower_rows[matrix.nrows - new_pivots + 1]
        end
    end

    matrix.npivots = matrix.nrows = matrix.size = new_pivots
    resize!(matrix.lower_rows, new_pivots)
    resize!(not_reduced_to_zero, new_pivots)
    true, any_zeroed, not_reduced_to_zero
end

function learn_interreduce_matrix_pivots!(
    graph::ComputationGraphF4,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A;
    reversed_rows::Bool=false
) where {C <: Coeff, A <: AbstractArithmetic}
    flag, _, not_reduced_to_zero =
        interreduce_matrix_pivots!(matrix, basis, arithmetic, reversed_rows=reversed_rows)
    !flag && return flag

    # Update the graph
    push!(
        graph.matrix_infos,
        (nup=matrix.nupper, nlow=matrix.nlower, ncols=matrix.ncolumns)
    )
    push!(graph.matrix_nonzeroed_rows, not_reduced_to_zero)
    push!(
        graph.matrix_upper_rows,
        (
            matrix.upper_to_coeffs[1:(matrix.nupper)],
            matrix.upper_to_mult[1:(matrix.nupper)]
        )
    )
    push!(graph.matrix_lower_rows, (Vector{Int}(), Vector{Int}()))
    true
end

function apply_interreduce_matrix_pivots!(
    graph::ComputationGraphF4,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A;
    reversed_rows::Bool=false
) where {C <: Coeff, A <: AbstractArithmetic}
    flag, any_zeroed, _ =
        interreduce_matrix_pivots!(matrix, basis, arithmetic, reversed_rows=reversed_rows)
    flag && !any_zeroed
end

# The number of blocks to split the matrix into
function nblocks_in_randomized(nrows::Int)
    floor(Int, sqrt(nrows / 3)) + 1
end

# Given a matrix of form 
#   A B
#   C D
# reduces the lower part CD with respect to the upper part AB.
# As a result, the matrix of the following form is produced:
#   A B
#   0 D' 
function randomized_reduce_matrix_lower_part!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A,
    rng
) where {C <: Coeff, A <: AbstractArithmetic}
    # Prepare the matrix
    pivots, row_idx_to_coeffs = absolute_index_pivots!(matrix)
    ncols = matrix.ncolumns
    nlow = matrix.nlower
    pivots = matrix.pivots
    row = zeros(C, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(C)

    nblocks = nblocks_in_randomized(nlow)
    rem = nlow % nblocks == 0 ? 0 : 1
    rowsperblock = div(nlow, nblocks) + rem
    @log level = -3 """
    Rows in the lower part: $nlow
    The bumber of blocks: $nblocks
    Rows per block: $rowsperblock"""
    nnz_counter = 0

    row = zeros(C, ncols)
    rows_multipliers = zeros(C, rowsperblock)

    for i in 1:nblocks
        nrowsupper = nlow > i * rowsperblock ? i * rowsperblock : nlow
        nrowstotal = nrowsupper - (i - 1) * rowsperblock
        nrowstotal == 0 && continue

        ctr = 0
        @inbounds while ctr < nrowstotal
            for j in 1:nrowstotal
                rows_multipliers[j] = mod_x(rand(rng, C), arithmetic)
            end
            row .= C(0)
            startcol = ncols

            for k in 1:nrowstotal
                rowidx = (i - 1) * rowsperblock + k
                nnz_column_indices = matrix.lower_rows[rowidx]
                nnz_coeffs = basis.coeffs[row_idx_to_coeffs[rowidx]]
                startcol = min(startcol, nnz_column_indices[1])

                for l in 1:length(nnz_coeffs)
                    colidx = nnz_column_indices[l]
                    row[colidx] =
                        mod_x(row[colidx] + rows_multipliers[k] * nnz_coeffs[l], arithmetic)
                end
            end

            # reduce it with known pivots from matrix.upper_rows
            # first nonzero in densecoeffs is at startcol position
            zeroed = reduce_dense_row_by_pivots_sparse!(
                new_column_indices,
                new_coeffs,
                row,
                matrix,
                basis,
                pivots,
                startcol,
                arithmetic,
                tmp_pos=-1
            )

            if zeroed
                ctr = nrowstotal
                break
            end
            nnz_counter += 1
            absolute_i = (i - 1) * rowsperblock + ctr + 1

            # matrix coeffs sparsely stores coefficients of new row
            matrix.coeffs[absolute_i] = new_coeffs
            # add new pivot at column index newrow[1]
            #  (which is the first nnz column of newrow)
            pivots[new_column_indices[1]] = new_column_indices
            # set ref to coefficient to matrix
            # guaranteed to be from lower part
            matrix.lower_to_coeffs[new_column_indices[1]] = absolute_i

            # normalize if needed
            normalize_sparse_row!(matrix.coeffs[absolute_i], arithmetic)

            ctr += 1
            new_column_indices, new_coeffs = new_empty_sparse_row(C)
        end
    end
    @log level = -3 "Nonzero rows: $nnz_counter (out of $nlow)"
    true
end

function reduce_matrix_lower_part_invariant_pivots!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    _, row_idx_to_coeffs = absolute_index_pivots!(matrix)
    ncols = matrix.ncolumns
    nlow = matrix.nlower
    pivots = matrix.pivots
    row = zeros(C, ncols)
    new_column_indices, new_coeffs = new_empty_sparse_row(C)
    @inbounds for i in 1:nlow
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = matrix.lower_rows[i]

        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref = matrix.coeffs[row_idx_to_coeffs[i]]
        load_sparse_row!(row, rowexps, cfsref)

        # reduce it with known pivots from matrix.upper_rows
        # first nonzero in densecoeffs is at startcol position
        first_nnz_column = rowexps[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            # first_nnz_column,
            # TODO: this is incorrect; for the counter-example see
            # https://github.com/sumiya11/Groebner.jl/issues/82
            1,
            arithmetic,
            tmp_pos=-1
        )

        # matrix coeffs sparsely stores coefficients of new row
        matrix.coeffs[i] = new_coeffs
        # add new pivot at column index newrow[1]
        #  (which is the first nnz column of newrow)
        matrix.lower_rows[i] = new_column_indices
        # set ref to coefficient to matrix
        # guaranteed to be from lower part
        matrix.lower_to_coeffs[i] = i

        new_column_indices, new_coeffs = new_empty_sparse_row(C)
    end
    matrix.npivots = matrix.nrows = matrix.size = matrix.nlower
end

function reduce_matrix_lower_part_any_nonzero!(
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    pivots, row_idx_to_coeffs = absolute_index_pivots!(matrix)
    row = zeros(C, matrix.ncolumns)
    new_column_indices, new_coeffs = new_empty_sparse_row(C)
    @inbounds for i in 1:(matrix.nlower)
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = matrix.lower_rows[i]
        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref = basis.coeffs[row_idx_to_coeffs[i]]
        load_sparse_row!(row, rowexps, cfsref)

        # reduce it with known pivots from matrix.upper_rows
        # first nonzero in densecoeffs is at startcol position
        first_nnz_column = rowexps[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            arithmetic,
            tmp_pos=-1
        )

        # if fully reduced
        zeroed && continue
        return false
    end
    true
end

function absolute_index_pivots!(matrix::MacaulayMatrix)
    pivots = Vector{Vector{ColumnLabel}}(undef, matrix.ncolumns)
    @inbounds for i in 1:(matrix.nupper)
        pivots[matrix.upper_rows[i][1]] = matrix.upper_rows[i]
    end
    absolute_lower_to_coeffs =
        Vector{ColumnLabel}(undef, max(matrix.ncolumns, matrix.nlower))
    @inbounds for i in 1:(matrix.nlower)
        absolute_lower_to_coeffs[matrix.lower_rows[i][1]] = matrix.lower_to_coeffs[i]
    end
    lower_to_coeffs = matrix.lower_to_coeffs
    matrix.lower_to_coeffs = absolute_lower_to_coeffs
    matrix.pivots = pivots
    pivots, lower_to_coeffs
end

function absolute_index_pivots_in_interreduction!(matrix::MacaulayMatrix, basis::Basis)
    resize!(matrix.lower_rows, matrix.ncolumns)
    resize!(matrix.upper_to_coeffs, matrix.ncolumns)
    resize!(matrix.upper_to_mult, matrix.ncolumns)
    resize!(matrix.lower_to_coeffs, matrix.ncolumns)
    resize!(matrix.lower_to_mult, matrix.ncolumns)
    resize!(matrix.coeffs, matrix.ncolumns)

    pivots = Vector{Vector{ColumnLabel}}(undef, matrix.ncolumns)
    @inbounds for i in 1:(matrix.nrows)
        pivots[matrix.upper_rows[i][1]] = matrix.upper_rows[i]
        matrix.lower_to_coeffs[matrix.upper_rows[i][1]] = i
        matrix.coeffs[i] = copy(basis.coeffs[matrix.upper_to_coeffs[i]])
    end

    matrix.pivots = pivots
    pivots
end

function new_empty_sparse_row(::Type{C}) where {C <: Coeff}
    Vector{ColumnLabel}(), Vector{C}()
end

# Normalize the row to have the leading coefficient equal to 1
function normalize_sparse_row!(
    row::Vector{T},
    arithmetic::A
) where {T <: CoeffFF, A <: AbstractArithmeticZp}
    @inbounds if isone(row[1])
        return row
    end
    @inbounds pinv = mod_x(invmod(row[1], divisor(arithmetic)), arithmetic)
    @inbounds row[1] = one(row[1])
    @inbounds for i in 2:length(row)
        row[i] = mod_x(row[i] * pinv, arithmetic)
    end
    row
end

# Normalize the row to have the leading coefficient equal to 1
function normalize_sparse_row!(
    row::Vector{T},
    arithmetic::A
) where {T <: CoeffQQ, A <: AbstractArithmeticQQ}
    @inbounds if isone(row[1])
        return row
    end
    @inbounds pinv = inv(row[1])
    @inbounds row[1] = one(row[1])
    @inbounds for i in 2:length(row)
        row[i] = row[i] * pinv
    end
    row
end

function reduce_dense_row_by_sparse_row!(
    row::Vector{T},
    indices::Vector{ColumnLabel},
    coeffs::Vector{T},
    arithmetic::A
) where {T <: CoeffFF, A <: AbstractArithmeticZp}
    @invariant isone(coeffs[1])
    @inbounds mul = divisor(arithmetic) - row[indices[1]]

    # On our benchmarks, usually,
    #   length(row) / length(indices)
    # roughly varies from 10 to 100
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = mod_x(row[idx] + mul * coeffs[j], arithmetic)
    end

    nothing
end

function reduce_dense_row_by_sparse_row!(
    row::Vector{T},
    indices::Vector{ColumnLabel},
    coeffs::Vector{T},
    arithmetic::A
) where {T <: CoeffQQ, A <: AbstractArithmeticQQ}
    @invariant isone(coeffs[1])

    @inbounds mul = -row[indices[1]]
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        row[idx] = row[idx] + mul * coeffs[j]
        # Base.GMP.MPQ.mul!(buf2, cfs[j], mul)
        # row[idx] = row[idx] + buf2
    end

    nothing
end

# Load coefficients from `coeffs` into dense `row` at `indices`
function load_sparse_row!(row::Vector{T}, indices, coeffs) where {T <: CoeffFF}
    @inbounds for i in 1:length(row)
        row[i] = T(0)
    end
    @inbounds for j in 1:length(indices)
        row[indices[j]] = coeffs[j]
    end
end

# Load coefficients from `coeffs` into dense `row` at `indices`
function load_sparse_row!(row::Vector{T}, indices, coeffs) where {T <: CoeffQQ}
    row .= T(0)
    @inbounds for j in 1:length(indices)
        row[indices[j]] = coeffs[j]
    end
end

# Traverses the dense `row` at positions `from..to` and extracts all nonzero
# entries to the given sparse row  
function extract_sparse_row!(indices, coeffs, row::Vector{T}, from::Int, to::Int) where {T}
    # NOTE: assumes the sparse row has the necessary capacity
    z = zero(T)
    j = 1
    @inbounds for i in from:to
        if row[i] != z
            indices[j] = i
            coeffs[j] = row[i]
            j += 1
        end
    end
    nothing
end

# function extract_sparse_row_2!(vals::Vector{UInt}, dense::Vector{UInt})
#     j = 1
#     z = Vec{4, UInt}((0, 0, 0, 0))
#     @inbounds for i in 1:8:length(dense)
#         v1 = SIMD.vload(Vec{4, UInt}, dense, i)
#         v2 = SIMD.vload(Vec{4, UInt}, dense, i + 4)
#         mask1 = !iszero(v1)
#         mask2 = !iszero(v2)
#         j1 = sum(mask1)
#         j2 = sum(mask2)
#         vstorec(v1, vals, j, mask1)
#         vstorec(v2, vals, j + j1, mask2)
#         j += j1 + j2
#     end
#     nothing
# end

record_active_reducer(active_reducers::Nothing, matrix, idx) = nothing
function record_active_reducer(active_reducers, matrix, idx)
    push!(active_reducers, (idx, matrix.upper_to_coeffs[idx], matrix.upper_to_mult[idx]))
    nothing
end

function reduce_dense_row_by_pivots_sparse!(
    new_column_indices::Vector{ColumnLabel},
    new_coeffs::Vector{C},
    row::Vector{C},
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    pivots::Vector{Vector{ColumnLabel}},
    start_column::Integer,
    arithmetic::A,
    active_reducers=nothing;
    tmp_pos::Integer=-1,
    exact_column_mapping::Bool=false
) where {C <: Coeff, A <: AbstractArithmetic}
    ncolumns = matrix.ncolumns
    nleft = matrix.nleft
    # The number of nonzeros in the reduced row
    nonzeros = 0
    # The index of the first nonzero elemen in the reduced row
    new_pivot = -1

    @inbounds for i in start_column:ncolumns
        # @inbounds for i in start_column:nleft
        # if the element is zero - no reduction is needed
        if iszero(row[i])
            continue
        end

        # if there is no pivot with the leading column index equal to i
        if !isassigned(pivots, i) || (tmp_pos != -1 && tmp_pos == i)
            if new_pivot == -1
                new_pivot = i
            end
            nonzeros += 1
            continue
        end
        # at this point, we have determined that row[i] is nonzero and that
        # pivots[i] is a valid reducer of row

        indices = pivots[i]
        if exact_column_mapping
            # if reducer is from new matrix pivots
            coeffs = matrix.coeffs[tmp_pos]
        elseif i <= nleft
            # if reducer is from the upper part of the matrix
            coeffs = basis.coeffs[matrix.upper_to_coeffs[i]]
            record_active_reducer(active_reducers, matrix, i)
        else
            # if reducer is from the lower part of the matrix
            coeffs = matrix.coeffs[matrix.lower_to_coeffs[i]]
        end

        reduce_dense_row_by_sparse_row!(row, indices, coeffs, arithmetic)
        @invariant iszero(row[i])
    end
    # form and return the resulting row in sparse format
    # all reduced to zero!
    if nonzeros == 0
        return true
    end
    resize!(new_column_indices, nonzeros)
    resize!(new_coeffs, nonzeros)
    extract_sparse_row!(
        new_column_indices,
        new_coeffs,
        row,
        convert(Int, start_column),
        ncolumns
    )
    false
end

function learn_reduce_matrix_lower_part!(
    graph::ComputationGraphF4,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    pivots, row_idx_to_coeffs = absolute_index_pivots!(matrix)
    ncols = matrix.ncolumns
    nlow = matrix.nlower
    row = zeros(C, ncols)

    not_reduced_to_zero = Vector{Int}()
    useful_reducers = Set{Tuple{Int, Int, MonomIdx}}()
    @log level = -6 "Low to coef" row_idx_to_coeffs matrix.lower_to_mult

    new_column_indices, new_coeffs = new_empty_sparse_row(C)
    @inbounds for i in 1:nlow
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = matrix.lower_rows[i]
        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref = basis.coeffs[row_idx_to_coeffs[i]]

        # we load coefficients into dense array
        # into rowexps indices
        load_sparse_row!(row, rowexps, cfsref)

        # reduce it with known pivots from matrix.upper_rows
        # first nonzero in densecoeffs is at startcol position
        reducers = Tuple{Int, Int, MonomIdx}[]
        first_nnz_column = rowexps[1]
        zeroed = reduce_dense_row_by_pivots_sparse!(
            new_column_indices,
            new_coeffs,
            row,
            matrix,
            basis,
            pivots,
            first_nnz_column,
            arithmetic,
            reducers,
            tmp_pos=-1
        )
        # if fully reduced
        zeroed && continue
        # NOTE: we are not adding reducers from lowrows!
        push!(not_reduced_to_zero, i)
        for rr in reducers
            push!(useful_reducers, rr)
        end

        @log level = -7 "Not zero" i row_idx_to_coeffs[i] matrix.lower_to_mult[i]

        # matrix coeffs sparsely stores coefficients of new row
        matrix.coeffs[i] = new_coeffs
        # add new pivot at column index newrow[1]
        #  (which is the first nnz column of newrow)
        pivots[new_column_indices[1]] = new_column_indices
        # set ref to coefficient to matrix
        # guaranteed to be from lower part
        matrix.lower_to_coeffs[new_column_indices[1]] = i

        # normalize if needed
        normalize_sparse_row!(matrix.coeffs[i], arithmetic)

        new_column_indices, new_coeffs = new_empty_sparse_row(C)
    end

    # NOTE: we sort reducers by their original position in the array of pivots.
    # That way, they are already sorted at the apply stage.
    useful_reducers_sorted = sort(collect(useful_reducers), by=reducer -> reducer[1])
    @log level = -7 "" useful_reducers_sorted
    push!(
        graph.matrix_infos,
        (nup=matrix.nupper, nlow=matrix.nlower, ncols=matrix.ncolumns)
    )
    push!(graph.matrix_nonzeroed_rows, not_reduced_to_zero)
    push!(
        graph.matrix_upper_rows,
        (map(f -> f[2], useful_reducers_sorted), map(f -> f[3], useful_reducers_sorted))
    )
    push!(
        graph.matrix_lower_rows,
        (row_idx_to_coeffs[not_reduced_to_zero], matrix.lower_to_mult[not_reduced_to_zero])
    )
    # push!(graph.matrix_upper_rows_buffers, map(f -> pivots[f[1]], useful_reducers_sorted))
    # push!(graph.matrix_lower_rows_buffers, matrix.lower_rows[not_reduced_to_zero])

    true
end

function apply_sparse_rref!(
    graph::ComputationGraphF4,
    ring,
    matrix::MacaulayMatrix{C},
    basis::Basis{C},
    rng
) where {C <: Coeff}
    ncols = matrix.ncolumns
    nlow = matrix.nlower

    arithmetic = select_arithmetic(matrix.coeffs, ring.ch)

    # move known matrix pivots,
    # no copy
    pivs, l2c_tmp = get_absolute_pivots!(matrix)
    @log level = -7 "absolute_pivots!" pivs l2c_tmp

    rowidx2coef = matrix.lower_to_coeffs
    matrix.lower_to_coeffs = l2c_tmp

    # unknown pivots,
    # we will modify them inplace when reducing by pivs
    upivs = matrix.lower_rows
    densecoeffs = zeros(C, ncols)

    @inbounds for i in 1:nlow
        # select next row to be reduced
        # npiv ~ exponents
        rowexps = upivs[i]

        # corresponding coefficients from basis
        # (no need to copy here)
        cfsref = basis.coeffs[rowidx2coef[i]]

        # we load coefficients into dense array
        # into rowexps indices
        load_sparse_row!(densecoeffs, rowexps, cfsref)

        # reduce it with known pivots from matrix.upper_rows
        # first nonzero in densecoeffs is at startcol position
        startcol = rowexps[1]
        zeroed, newrow, newcfs = reduce_dense_row_by_pivots_sparse!(
            densecoeffs,
            matrix,
            basis,
            pivs,
            ColumnLabel(startcol),
            ColumnLabel(-1),
            arithmetic
        )
        # if fully reduced
        if zeroed
            # then something has gone wrong in tracing 
            return false
        end

        # matrix coeffs sparsely stores coefficients of new row
        matrix.coeffs[i] = newcfs
        # add new pivot at column index newrow[1]
        #  (which is the first nnz column of newrow)
        pivs[newrow[1]] = newrow
        # set ref to coefficient to matrix
        # guaranteed to be from lower part
        matrix.lower_to_coeffs[newrow[1]] = i

        # normalize if needed
        normalize_sparse_row!(matrix.coeffs[i], arithmetic)
    end

    interreduce_lower_part!(matrix, basis, pivs, arithmetic)
    true
end

function column_to_monom_mapping!(matrix::MacaulayMatrix, symbol_ht::MonomialHashtable)

    # monoms from symbolic table represent one column in the matrix
    hdata = symbol_ht.hashdata
    load = symbol_ht.load

    column_to_monom = Vector{MonomIdx}(undef, load - 1)
    j = 1
    # number of pivotal cols
    k = 0
    @inbounds for i in (symbol_ht.offset):load
        # column to hash index
        column_to_monom[j] = i
        j += 1
        # meaning the column is pivoted
        if hdata[i].idx == 2
            k += 1
        end
    end

    sort_columns_by_labels!(column_to_monom, symbol_ht)

    matrix.nleft = k  # CHECK!
    # -1 as long as hashtable load is always 1 more than actual
    matrix.nright = load - matrix.nleft - 1

    # store the other direction of mapping,
    # hash -> column
    @inbounds for k in 1:length(column_to_monom)
        hv = hdata[column_to_monom[k]]
        hdata[column_to_monom[k]] = Hashvalue(k, hv.hash, hv.divmask, hv.deg)
    end

    @inbounds for k in 1:(matrix.nupper)
        row = matrix.upper_rows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    @inbounds for k in 1:(matrix.nlower)
        row = matrix.lower_rows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    matrix.ncolumns = matrix.nleft + matrix.nright

    @assert matrix.nleft + matrix.nright == symbol_ht.load - 1 == matrix.ncolumns
    @assert matrix.nlower + matrix.nupper == matrix.nrows

    matrix.column_to_monom = column_to_monom
end

function convert_rows_to_basis_elements!(
    matrix::MacaulayMatrix,
    basis::Basis,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable
)

    # we mutate basis array directly by adding new elements

    resize_basis_if_needed!(basis, matrix.npivots)
    rows = matrix.lower_rows
    crs = basis.nprocessed

    @inbounds for i in 1:(matrix.npivots)
        colidx = rows[i][1]
        insert_in_basis_hash_table_pivots(rows[i], ht, symbol_ht, matrix.column_to_monom)
        basis.coeffs[crs + i] = matrix.coeffs[matrix.lower_to_coeffs[colidx]]
        basis.monoms[crs + i] = matrix.lower_rows[i]
    end

    basis.nfilled += matrix.npivots
end

function convert_rows_to_basis_elements_nf!(
    matrix::MacaulayMatrix,
    basis::Basis,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable
)
    resize_basis_if_needed!(basis, matrix.npivots)

    @inbounds for i in 1:(matrix.npivots)
        basis.nprocessed += 1
        basis.nnonredundant += 1
        basis.nonredundant[basis.nnonredundant] = basis.nprocessed
        if isassigned(matrix.coeffs, i)
            row = matrix.lower_rows[i]
            insert_in_basis_hash_table_pivots(row, ht, symbol_ht, matrix.column_to_monom)
            basis.coeffs[basis.nprocessed] = matrix.coeffs[i]
            basis.monoms[basis.nprocessed] = row
        else
            empty!(basis.coeffs[basis.nprocessed])
            empty!(basis.monoms[basis.nprocessed])
        end
    end

    nothing
end

function insert_in_basis_hash_table_pivots(
    row::Vector{ColumnLabel},
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    column_to_monom::Vector{MonomIdx}
) where {M}
    resize_hashtable_if_needed!(ht, length(row))

    sdata = symbol_ht.hashdata
    sexps = symbol_ht.monoms

    mod = MonomHash(ht.size - 1)
    bdata = ht.hashdata
    bexps = ht.monoms
    bhash = ht.hashtable

    l = 1
    @label Letsgo
    @inbounds while l <= length(row)
        hidx = column_to_monom[row[l]]

        # symbolic hash
        h = sdata[hidx].hash

        lastidx = ht.load + 1
        bexps[lastidx] = sexps[hidx]
        e = bexps[lastidx]

        k = h
        i = MonomHash(1)
        @inbounds while i <= ht.size
            k = next_lookup_index(h, i, mod)
            hm = bhash[k]

            iszero(hm) && break

            if ishashcollision(ht, hm, e, h)
                i += MonomHash(1)
                continue
            end

            row[l] = hm
            l += 1
            @goto Letsgo
        end

        bhash[k] = pos = lastidx
        row[l] = pos
        l += 1

        bdata[pos] = Hashvalue(sdata[hidx].idx, h, sdata[hidx].divmask, sdata[hidx].deg)

        ht.load += 1
    end

    nothing
end

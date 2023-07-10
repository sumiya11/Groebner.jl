
function reduction_apply!(
    graph::ComputationGraphF4,
    ring,
    basis,
    matrix,
    ht,
    rng,
    symbol_ht
)
    column_to_monom_mapping!(matrix, symbol_ht)

    sort_matrix_upper_rows_decreasing!(matrix) # for pivots,  AB part
    sort_matrix_lower_rows_increasing!(matrix) # for reduced, CD part

    @log "Apply: after mapping and sorting" matrix
    #=
    ┌ Info: Apply: after mapping and sorting
    └   matrix = Main.Groebner.MacaulayMatrix{UInt64}(
        Vector{Int32}[[3, 2]], 
        Vector{Int32}[[3, 1]], 
        Int32[2, 4, 3], 
        Vector{UInt64}[], 
        2, 0, 2, 3, 1, 1, 1, 2, [1], 
    [2], Int32[], Int32[])
    =#
    flag = linear_algebra!(graph, ring, matrix, basis, :apply, rng)
    if !flag
        return false
    end
    @log "After linear algebra" matrix

    convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht)
    true
end

function symbolic_preprocessing!(
    graph::ComputationGraphF4,
    iters::Int,
    basis,
    matrix,
    hashtable,
    symbol_ht
)
    lowrows, lowmults = graph.matrix_lower_rows[iters]
    uprows, upmults = graph.matrix_upper_rows[iters]
    matrix_info = graph.matrix_infos[iters]
    nonzeroed_rows = graph.matrix_nonzeroed_rows[iters]
    nlow = length(nonzeroed_rows) # matrix_info.nlow # length(nonzeroed_rows)
    nup = length(uprows) # matrix_info.nup # length(uprows)
    # TODO: check how many columns are there
    # ncols = matrix_info.ncols

    matrix.uprows = Vector{Vector{ColumnIdx}}(undef, nup)
    matrix.lowrows = Vector{Vector{ColumnIdx}}(undef, nlow)
    matrix.low2coef = Vector{Int}(undef, nlow)
    matrix.up2coef = Vector{Int}(undef, nup)

    check_enlarge_hashtable!(symbol_ht, nlow + nup + 2)
    for i in 1:nlow
        mult_idx = lowmults[i]
        poly_idx = lowrows[i]

        h = hashtable.hashdata[mult_idx].hash
        etmp = hashtable.monoms[mult_idx]
        rpoly = basis.monoms[poly_idx]

        @log "Preserved" poly_idx hashtable.monoms[rpoly]

        @log "Low row" etmp

        # vidx = insert_in_hash_table!(symbol_ht, etmp)

        matrix.lowrows[i] =
            multiplied_poly_to_matrix_row!(symbol_ht, hashtable, h, etmp, rpoly)
        symbol_ht.hashdata[matrix.lowrows[i][1]].idx = 2  # TODO!!

        matrix.low2coef[i] = poly_idx

        @log "low row x2: " matrix.lowrows[i] symbol_ht.monoms[matrix.lowrows[i]]
    end

    for i in 1:nup
        mult_idx = upmults[i]
        poly_idx = uprows[i]

        h = hashtable.hashdata[mult_idx].hash
        etmp = hashtable.monoms[mult_idx]
        rpoly = basis.monoms[poly_idx]

        @log "Up row" etmp

        # vidx = insert_in_hash_table!(symbol_ht, etmp)

        matrix.uprows[i] =
            multiplied_poly_to_matrix_row!(symbol_ht, hashtable, h, etmp, rpoly)
        symbol_ht.hashdata[matrix.uprows[i][1]].idx = 2  # TODO!!

        matrix.up2coef[i] = poly_idx

        @log "up row x2: " matrix.uprows[i] symbol_ht.monoms[matrix.uprows[i]]
    end

    i = MonomIdx(symbol_ht.offset)
    @inbounds while i <= symbol_ht.load
        # not a reducer
        if iszero(symbol_ht.hashdata[i].idx)
            symbol_ht.hashdata[i].idx = 1
        end
        i += MonomIdx(1)
    end

    @log "Symbol ht:" symbol_ht

    matrix.nrows = nlow + nup
    matrix.nlow = nlow
    matrix.nup = nup
    matrix.size = matrix.nrows
end

function reducegb_f4_apply!(
    graph,
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    hashtable::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M}
) where {M}
    @log level = 100000 "Entering apply autoreduction" basis
    #=
    ┌ LogLevel(100000): Entering apply autoreduction
    │   basis = Main.Groebner.Basis{UInt64}(
    Vector{Int32}[[3], [2]], 
    Vector{UInt64}[[0x0000000000000001], 
    [0x0000000000000001]], 2, 2, 2, Bool[0, 0], 
    [1, 2], UInt32[0x00010000, 0x00000001], 2)    =#

    # TODO: check that reduced flag is set

    lowrows, lowmults = graph.matrix_lower_rows[end]
    uprows, upmults = graph.matrix_upper_rows[end]
    matrix_info = graph.matrix_infos[end]
    nonzeroed_rows = graph.matrix_nonzeroed_rows[end]
    nlow = length(lowrows) # matrix_info.nlow
    nup = length(uprows) # matrix_info.nup # 
    # ncols = matrix_info.ncols
    nonred = graph.nonredundant_indices_before_reduce

    matrix.uprows = Vector{Vector{ColumnIdx}}(undef, nup)
    matrix.lowrows = Vector{Vector{ColumnIdx}}(undef, nlow)
    matrix.low2coef = Vector{Int}(undef, nlow)
    matrix.up2coef = Vector{Int}(undef, nup)

    matrix.nrows = 0
    matrix.ncols = 0
    matrix.nleft = 0
    matrix.nright = 0
    matrix.nup = nup
    matrix.nlow = nlow

    etmp = construct_const_monom(M, hashtable.nvars)
    # etmp is now set to zero, and has zero hash

    check_enlarge_hashtable!(symbol_ht, nlow + nup + 2)
    # @inbounds for i in 1:length(nonred) #
    #     matrix.nrows += 1
    #     matrix.uprows[matrix.nrows] = multiplied_poly_to_matrix_row!(
    #         symbol_ht,
    #         hashtable,
    #         MonomHash(0),
    #         etmp,
    #         basis.monoms[nonred[i]]
    #     )

    #     matrix.up2coef[matrix.nrows] = nonred[i]
    #     # matrix.up2mult[matrix.nrows] = insert_in_hash_table!(hashtable, etmp)
    #     # set lead index as 1
    #     symbol_ht.hashdata[matrix.uprows[matrix.nrows][1]].idx = 1
    # end

    # needed for correct column count in symbol hashtable
    matrix.ncols = matrix.nrows
    matrix.nup = matrix.nrows

    @log "Before autoreduce apply" basis uprows upmults lowmults matrix_info
#=
   ┌ Info: Before autoreduce apply
│   basis = Main.Groebner.Basis{UInt64}(Vector{Int32}[[4, 5], [2, 3], [6, 3], #undef], Vector{UInt64}[[0x0000000000000001, 0x0000000000000001], [0x0000000000000001, 0x0000000000000001], [0x0000000000000001, 0x000000007ffffffe], #undef], 4, 3, 3, Bool[0, 0, 0, 0], [1, 2, 3, 3], UInt32[0x00000001, 0x00010001, 0x00030000, 0x00000000], 3)
│   uprows =
│    5-element Vector{Int64}:
│     2
│     3
│     1
│     3
│     3
│   upmults =
│    5-element Vector{Int32}:
│              3
│              3
│              3
│            572
│     1714710608
│   lowmults = Int32[]
└   matrix_info = (nup = 3, nlow = 0, ncols = 5)
=#
    for i in 1:nup
        mult_idx = upmults[i]
        poly_idx = uprows[i]

        h = hashtable.hashdata[mult_idx].hash
        etmp = hashtable.monoms[mult_idx]
        rpoly = basis.monoms[poly_idx]

        @log "Up row" etmp

        # Not too good
        # iszero(h) && continue

        # vidx = insert_in_hash_table!(symbol_ht, etmp)
        matrix.uprows[i] =
            multiplied_poly_to_matrix_row!(symbol_ht, hashtable, h, etmp, rpoly)
        # symbol_ht.hashdata[matrix.uprows[i][1]].idx = 2  # TODO!!

        matrix.up2coef[i] = poly_idx

        # @log "up row x2: " matrix.uprows[matrix.nrows] symbol_ht.monoms[matrix.uprows[matrix.nrows]]
    end

    # set all pivots to unknown
    @inbounds for i in (symbol_ht.offset):(symbol_ht.load)
        symbol_ht.hashdata[i].idx = 1
    end
    # matrix.ncols = ncols

    matrix.nrows = nlow + nup
    matrix.nlow = nlow
    matrix.nup = nup
    matrix.size = matrix.nrows

    column_to_monom_mapping!(matrix, symbol_ht)
    matrix.ncols = matrix.nleft + matrix.nright

    sort_matrix_upper_rows_decreasing!(matrix)
    
    @log "In autoreduction apply" basis matrix

#=
┌ Info: In autoreduction apply
│   basis = Main.Groebner.Basis{UInt64}(
Vector{Int32}[[4, 5], [2, 3], [6, 3], 
#undef], 
Vector{UInt64}[[0x0000000000000001, 0x0000000000000001], [0x0000000000000001, 0x0000000000000001], [0x0000000000000001, 0x000000007ffffffe], #undef], 
4, 3, 3, Bool[0, 0, 0, 0], [1, 2, 3, 3], UInt32[0x00000001, 0x00010001, 0x00030000, 0x00000000], 3)
└   matrix = Main.Groebner.MacaulayMatrix{UInt64}(
    Vector{Int32}[[1, 5], [2, 5], [3, 4]], 
    Vector{Int32}[], 
    Int32[2, 4, 5, 6, 3], 
    Vector{UInt64}[], 
    3, 0, 3, 5, 3, 0, 0, 5, 
    [2, 3, 1], 
    Int64[], Int32[], Int32[])
=#

    flag = exact_sparse_rref_interreduce_apply!(ring, matrix, basis)
    if !flag
        return false
    end
    convert_rows_to_basis_elements!(matrix, basis, hashtable, symbol_ht)

    basis.ntotal = matrix.npivots + basis.nprocessed
    basis.nprocessed = matrix.npivots

    @log level = 100000 "Apply autoreduction: after linalg" basis matrix
    #=
    ┌ LogLevel(100000): Apply autoreduction: after linalg
    │   basis = Main.Groebner.Basis{UInt64}(
    Vector{Int32}[[2], [2], [2], #undef], 
    Vector{UInt64}[[0x0000000000000001], [0x0000000000000001], 
    [0x0000000000000001], #undef], 
    4, 1, 3, 
    Bool[0, 0, 0, 0], 
    [1, 2, 563204, 3], 
    UInt32[0x00000001, 0x00000001, 0x9cfd5eb0, 0x00000143], 
    2)    
    │   matrix = Main.Groebner.MacaulayMatrix{UInt64}(
    Vector{Int32}[[1]], Vector{Int32}[[2]], 
    Int32[2], 
    Vector{UInt64}[[0x0000000000000001]], 
    1, 1, 1, 1, 1, 0, 0, 1, 
    [2], [1], 
    Int32[178984752], Int32[2085144128])
    └ @ Main.Groebner c:\data\projects\gbgb\Groebner.jl\src\utils\logging
    =#
    # we may have added some multiples of reduced basis polynomials
    # from the matrix, so get rid of them
    output_nonredundant = graph.output_nonredundant_indices
    # basis.nonredundant = output_nonredundant
    for i in 1:length(output_nonredundant)
        basis.nonredundant[i] = output_nonredundant[i]
        basis.divmasks[i] =
            hashtable.hashdata[basis.monoms[basis.nonredundant[i]][1]].divmask
    end
    basis.ndivmasks = length(output_nonredundant)
    true
end

function standardize_basis_apply!(graph)
    basis = graph.gb_basis
    buf = graph.buf_basis
    basis.size = basis.nprocessed = basis.ntotal = basis.ndivmasks
    output_nonredundant = graph.output_nonredundant_indices
    output_sort = graph.output_sort_indices
    @inbounds for i in 1:basis.nprocessed
        basis.coeffs[i] = buf.coeffs[output_nonredundant[output_sort[i]]]
    end
    buf.nprocessed = buf.ndivmasks = 0
    buf.ntotal = graph.input_basis.ntotal
    normalize_basis!(graph.ring, basis)
end

function standardize_basis_learn!(graph, ring, basis, ht, ord)
    @inbounds for i in 1:(basis.ndivmasks)
        idx = basis.nonredundant[i]
        basis.nonredundant[i] = i
        basis.isredundant[i] = false
        basis.coeffs[i] = basis.coeffs[idx]
        basis.monoms[i] = basis.monoms[idx]
    end
    basis.size = basis.nprocessed = basis.ntotal = basis.ndivmasks
    resize!(basis.coeffs, basis.nprocessed)
    resize!(basis.monoms, basis.nprocessed)
    resize!(basis.divmasks, basis.nprocessed)
    resize!(basis.nonredundant, basis.nprocessed)
    resize!(basis.isredundant, basis.nprocessed)
    perm = sort_polys_by_lead_increasing!(basis, ht, ord=ord)
    graph.output_sort_indices = perm
    normalize_basis!(ring, basis)
end

function f4_apply!(graph, ring, basis::Basis{C}, params) where {M <: Monom, C <: Coeff}
    @invariant basis_well_formed(:input_f4_apply!, ring, basis, graph.hashtable)
    @assert params.reduced === true

    normalize_basis!(ring, basis)

    iters_total = length(graph.matrix_infos) - 1
    iters = 0
    hashtable = graph.hashtable

    matrix = initialize_matrix(ring, C)
    @log "Applying modulo $(ring.ch)"

    update_basis!(basis, hashtable)

    @log level = 100 "Input basis:" basis
    #=
   ┌ LogLevel(100): Input basis:
└   basis = Main.Groebner.Basis{UInt64}(
Vector{Int32}[[4, 5], [2, 3]], 
Vector{UInt64}[
[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x0000000000000001]], 
2, 2, 2, 
Bool[0, 0], [1, 2], 
UInt32[0x00000001, 0x00010001], 
2) 
=#

    while iters < iters_total
        iters += 1
        @log "F4 Apply iteration $iters"

        symbol_ht = initialize_secondary_hashtable(hashtable)

        symbolic_preprocessing!(graph, iters, basis, matrix, hashtable, symbol_ht)
        @log "After symbolic preprocessing:" matrix
        #=
        
        =#

        flag = reduction_apply!(graph, ring, basis, matrix, hashtable, params.rng, symbol_ht)
        if !flag
            # Unlucky cancellation of basis coefficients happened
            return false
        end
        #=
        ┌ Info: After reduction_apply:
        │   matrix = Main.Groebner.MacaulayMatrix{UInt64}(
        Vector{Int32}
        [[1, 3, 4], [2, 4, 5]], 
        Vector{Int32}[[10, 7, 11]], 
        Int32[2, 3, 5, 4, 6], 
        Vector{UInt64}[[0x0000000000000001, 0x0000000000000001, 0x0000000000000001]], 
        1, 1, 1, 5, 2, 1, 2, 3, 
        [1, 1], 
        [2, -1, 1, 32766, -321781503], Int32[], Int32[])

        └   basis = Main.Groebner.Basis{UInt64}(
        Vector{Int32}[
        [2, 3, 4], [5, 6, 7], 
        [8, 9], [10, 7, 11], #undef, #undef], 
        Vector{UInt64}[[0x0000000000000001, 0x0000000000000001, 0x0000000000000001], [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], [0x0000000000000001, 0x0000000000000010], [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], #undef, #undef], 
        6, 3, 4, 
        Bool[0, 0, 0, 0, 0, 0], 
        [1, 2, 3, 4, 4, 1438045129312], 
        UInt32[0x00000001, 0x00000401, 0x00100401, 0x0000014e, 0xc630f010, 0x00007ffe], 
        3)
        =#

        @log "After reduction_apply:" matrix basis

        update_basis!(basis, hashtable)

        @log "After update apply" basis

        #=
 ┌ Info: After update apply
└   basis = Main.Groebner.Basis{UInt64}(
Vector{Int32}[
    [4, 5], [2, 3], [6, 3], #undef], 
Vector{UInt64}
[[0x0000000000000001, 0x0000000000000001],
[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x000000007ffffffe], 
#undef], 
4, 3, 3, 
Bool[0, 0, 0, 0], 
[1, 2, 3, 3], 
UInt32[0x00000001, 0x00010001, 0x00030000, 0x00000000], 
3)

┌ Info: After update apply
└   basis = Main.Groebner.Basis{UInt64}(
Vector{Int32}[[4, 5], [2, 3], [6, 3], 
#undef], 
Vector{UInt64}[[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x000000007ffffffe], #undef], 
4, 3, 3, Bool[0, 0, 0, 0], [1, 2, 3, 3], 
UInt32[0x00000001, 0x00010001, 0x00030000, 0x00000000], 
3)
        =#

        matrix = initialize_matrix(ring, C)
    end

    @log "Before reduction" basis


    # filter_redundant!(basis)

    if params.reduced
        # TODO: drop redundant!!
        @log "Autoreducing the final basis.."
        symbol_ht = initialize_secondary_hashtable(hashtable)
        flag = reducegb_f4_apply!(graph, ring, basis, matrix, hashtable, symbol_ht)
        if !flag
            return false
        end
        @log "Autoreduced!"
    end

    
    standardize_basis_apply!(graph)
    basis = graph.gb_basis

    @log "After apply standardization" basis
    #=
    ┌ Info: After apply standardization
    └   basis = Main.Groebner.Basis{UInt64}(
    Vector{Int32}[[2]], 
    Vector{UInt64}[[0x0000000000000001]], 
    1, 1, 1, 
    Bool[0], [1], 
    UInt32[0x00000001], 
    1) 
    =#

    @invariant basis_well_formed(:output_f4_apply!, ring, basis, hashtable)

    true
end

function update!(graph, pairset, basis, hashtable, update_ht)
    update!(pairset, basis, hashtable, update_ht)
end

function symbolic_preprocessing!(graph, basis, matrix, hashtable, symbol_ht)
    symbolic_preprocessing!(basis, matrix, hashtable, symbol_ht)
    @log "Symbol ht:" symbol_ht
end

function reduction_learn!(graph, ring, basis, matrix, hashtable, symbol_ht, linalg, rng)
    column_to_monom_mapping!(matrix, symbol_ht)

    sort_matrix_upper_rows_decreasing!(matrix) # for pivots,  AB part
    sort_matrix_lower_rows_increasing!(matrix) # for reduced, CD part

    @log "Learn: after mapping and sorting" matrix
    #=

    ┌ Info: Learn: after mapping and sorting
    └   matrix = Main.Groebner.MacaulayMatrix{UInt64}(
    Vector{Int32}[[1, 2]], 
    Vector{Int32}[[1, 3]], 
    Int32[2, 3, 4], 
    Vector{UInt64}[], 
    2, 0, 2, 3, 1, 1, 1, 2, 
    [1, 2, 559108, 2, 2, 1670339522160, 1, 0],
     [2, 2], 
     Int32[6, 388, 1, 0, 563205, 0, 1, 0],
    Int32[7, 388])    
    =#

    linear_algebra!(graph, ring, matrix, basis, :learn, rng)

    @log "After linear algebra" matrix

    convert_rows_to_basis_elements!(matrix, basis, hashtable, symbol_ht)
end

function f4_reducegb_learn!(
    graph,
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M}
) where {M}
    @log level = 100000 "Entering learn autoreduction" basis

    #=
  
    =#

    etmp = construct_const_monom(M, ht.nvars)
    # etmp is now set to zero, and has zero hash

    reinitialize_matrix!(matrix, basis.ndivmasks)
    uprows = matrix.uprows

    @log "Before autoreduce learn" basis

    # add all non redundant elements from basis
    # as matrix upper rows
    @inbounds for i in 1:(basis.ndivmasks) #
        matrix.nrows += 1
        uprows[matrix.nrows] = multiplied_poly_to_matrix_row!(
            symbol_ht,
            ht,
            MonomHash(0),
            etmp,
            basis.monoms[basis.nonredundant[i]]
        )

        matrix.up2coef[matrix.nrows] = basis.nonredundant[i]
        matrix.up2mult[matrix.nrows] = insert_in_hash_table!(ht, etmp)
        # set lead index as 1
        symbol_ht.hashdata[uprows[matrix.nrows][1]].idx = 1
    end
    # TODO
    graph.nonredundant_indices_before_reduce = basis.nonredundant[1:(basis.ndivmasks)]

    # needed for correct column count in symbol hashtable
    matrix.ncols = matrix.nrows
    matrix.nup = matrix.nrows

    symbolic_preprocessing!(basis, matrix, ht, symbol_ht)
    # set all pivots to unknown
    @inbounds for i in (symbol_ht.offset):(symbol_ht.load)
        symbol_ht.hashdata[i].idx = 1
    end

    column_to_monom_mapping!(matrix, symbol_ht)
    matrix.ncols = matrix.nleft + matrix.nright

    sort_matrix_upper_rows_decreasing!(matrix)

    @log "In autoreduction learn" basis matrix
#=
┌ Info: In autoreduction learn
│   basis = Main.Groebner.Basis{UInt64}(Vector{Int32}
[[4, 5], [2, 3], [6, 3], #undef], 
Vector{UInt64}[[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x000000007ffffffe], #undef], 
4, 3, 3, Bool[0, 0, 0, 0], [1, 2, 3, 19], 
UInt32[0x00000001, 0x00010001, 0x00030000, 0x00000000], 
3)
└   matrix = Main.Groebner.MacaulayMatrix{UInt64}(
Vector{Int32}[[1, 5], [2, 5], [3, 4]], 
Vector{Int32}[#undef, #undef, #undef, #undef, #undef, #undef], 
Int32[4, 6, 2, 3, 5], Vector{UInt64}[], 
3, 0, 3, 5, 3, 0, 0, 5, 
[2, 3, 1, -1, 0, 7488, 7, 8, 9, 10, 12, 0], 
[48, 3, 563204, 3, 3, 2458435584008], 
Int32[3, 3, 3, 0, 563205, 0, 982866448, 572, 982865936, 572, 982865808, 572], 
Int32[4, 0, 2112605776, 32767, 563200, 0])   
=#

    exact_sparse_rref_interreduce_learn!(graph, ring, matrix, basis)

    convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht)

    basis.ntotal = matrix.npivots + basis.nprocessed
    basis.nprocessed = matrix.npivots

    @log level = 100000 "Learn autoreduction: after linalg" basis matrix

    #=
│   basis = Main.Groebner.Basis{UInt64}(
Vector{Int32}[[4, 5], [2, 3], [6, 3], 
[2, 3], [6, 3], [4, 5], 
#undef, #undef], 
Vector{UInt64}[[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x000000007ffffffe], 
[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x000000007ffffffe], 
[0x0000000000000001, 0x0000000000000001], 
#undef, #undef], 
8, 3, 6, 
Bool[0, 0, 0, 0, 0, 0, 0, 0], 
[1, 2, 3, 19, 22, 23, 38, 39], 
UInt32[0x00000001, 0x00010001, 0x00030000, 0x00000000, 0x00089805, 0x00000000, 0x0000001e, 
0x00000000], 3)
    =#

    # we may have added some multiples of reduced basis polynomials
    # from the matrix, so get rid of them
    k = 0
    i = 1
    @label Letsgo
    @inbounds while i <= basis.nprocessed
        @inbounds for j in 1:k
            if is_monom_divisible(
                basis.monoms[basis.ntotal - i + 1][1],
                basis.monoms[basis.nonredundant[j]][1],
                ht
            )
                i += 1
                @goto Letsgo
            end
        end
        k += 1
        basis.nonredundant[k] = basis.ntotal - i + 1
        basis.divmasks[k] = ht.hashdata[basis.monoms[basis.nonredundant[k]][1]].divmask
        i += 1
    end
    basis.ndivmasks = k

    graph.output_nonredundant_indices = copy(basis.nonredundant[1:k])
end

function f4_learn!(
    graph,
    ring::PolyRing,
    basis::Basis{C},
    pairset::Pairset,
    hashtable::MonomialHashtable{M},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    # @invariant hashtable_well_formed(:input_f4!, ring, hashtable)
    @invariant basis_well_formed(:input_f4_learn!, ring, basis, hashtable)
    # @invariant pairset_well_formed(:input_f4!, pairset, basis, ht)

    @assert basis == graph.gb_basis
    @assert params.reduced === true

    @log level = 5 "Entering F4 Learn phase."
    # TODO: decide on the number field arithmetic implementation
    normalize_basis!(ring, basis)

    matrix = initialize_matrix(ring, C)

    # initialize hash tables for update and symbolic preprocessing steps
    update_ht = initialize_secondary_hashtable(hashtable)
    symbol_ht = initialize_secondary_hashtable(hashtable)

    # add the first batch of critical pairs to the pairset
    @log level = 6 "Processing initial polynomials, generating first critical pairs"
    pairset_size = update!(graph, pairset, basis, hashtable, update_ht)
    @log level = 6 "Out of $(basis.ntotal) polynomials, $(basis.nprocessed) are non-redundant"
    @log level = 6 "Generated $(pairset.load) critical pairs"

    # TODO: the order of input polynomials is not deterministic!!
    # We need to remember the permutation of the input polynomials
    @log level = 100 "Input basis:" basis

    #=
└   basis = Main.Groebner.Basis{UInt64}(
   basis = Main.Groebner.Basis{UInt64}(
Vector{Int32}[[4, 5], [2, 3]], 
Vector{UInt64}[[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x0000000000000001]], 
2, 2, 2, Bool[0, 0], [1, 2], UInt32[0x00000001, 0x00010001], 2) 
    =#
    i = 0
    # While there are pairs to be reduced
    while !isempty(pairset)
        i += 1
        @log "F4: iteration $i"
        @log "F4: available $(pairset.load) pairs"

        # selects pairs for reduction from pairset following normal strategy
        # (minimal lcm degrees are selected),
        # and puts these into the matrix rows
        select_normal!(
            pairset,
            basis,
            matrix,
            hashtable,
            symbol_ht,
            maxpairs=params.maxpairs
        )
        # Color with [F4]

        symbolic_preprocessing!(graph, basis, matrix, hashtable, symbol_ht)
        @log "Formed a matrix of size X, DISPLAY_MATRIX"
        @log "After symbolic preprocessing:" matrix

#=
┌ Info: After symbolic preprocessing:
└   matrix = Main.Groebner.MacaulayMatrix{UInt64}(
Vector{Int32}[[2, 3]], 
Vector{Int32}[[2, 4]], Int32[], Vector{UInt64}[], 
2, 0, 2, 3, 1, 1, 0, 0, 
[1, 0, 0, 0, 0, 0, 0, 0], [2, 0], 
Int32[5, 0, 0, 0, 0, 0, 0, 0], 
Int32[3, 0])
=#

#=
┌ Info: After symbolic preprocessing:
└   matrix = Main.Groebner.MacaulayMatrix{UInt64}(
Vector{Int32}[[2, 3], [4, 3]], 
Vector{Int32}[[2, 4]], 
Int32[], Vector{UInt64}[], 
3, 0, 3, 3, 2, 1, 0, 0, 
[2, 1, 563204, 3, 3, 140735303642192, 2456993986752, 2456993987296], 
[3, 3], 
Int32[5, 3, 1714506584, 572, 1714729136, 572, 565, 1], 
Int32[4, 32767])
=#

        # reduces polys and obtains new potential basis elements
        reduction_learn!(
            graph,
            ring,
            basis,
            matrix,
            hashtable,
            symbol_ht,
            params.linalg,
            params.rng
        )
        #=
┌ Info: After reduction_learn:
│   matrix = Main.Groebner.MacaulayMatrix{UInt64}(
Vector{Int32}[[1, 2]], Vector{Int32}[[6, 3]], 
Int32[2, 3, 4], 
Vector{UInt64}[[0x0000000000000001, 0x000000007ffffffe]], 
1, 1, 1, 3, 1, 1, 1, 2, 
[1, 0, 0, 0, 0, 0, 0, 0], 
[2, 1, 0], 
Int32[5, 0, 0, 0, 0, 0, 0, 0], 
Int32[3, 0])        
└   basis = Main.Groebner.Basis{UInt64}(
Vector{Int32}[[4, 5], [2, 3], [6, 3], #undef], 
Vector{UInt64}[[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x000000007ffffffe], #undef], 
4, 2, 3, 
Bool[0, 0, 0, 0], [1, 2, 18, 19], 
UInt32[0x00000001, 0x00010001, 0x0000001e, 0x00000000], 
2)
        =#

        # TODO: insert monoms from symbol_ht into the main hashtable
        @log "After reduction_learn:" matrix basis

        # update the current basis with polynomials produced from reduction,
        # does not copy,
        # checks for redundancy
        pairset_size = update!(graph, pairset, basis, hashtable, update_ht)
        # TODO: move the above

        @log "After update learn" basis
        #=
┌ Info: After update learn
└   basis = Main.Groebner.Basis{UInt64}(
Vector{Int32}[[4, 5], [2, 3], [6, 3], #undef], 
Vector{UInt64}[[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x000000007ffffffe], #undef], 
4, 3, 3, Bool[0, 0, 0, 0], [1, 2, 3, 19], 
UInt32[0x00000001, 0x00010001, 0x00030000, 0x00000000], 
3)

┌ Info: After update learn
└   basis = Main.Groebner.Basis{UInt64}(
Vector{Int32}[[4, 5], [2, 3], [6, 3], #undef], 
Vector{UInt64}[[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x0000000000000001], 
[0x0000000000000001, 0x000000007ffffffe], #undef],
 4, 3, 3, 
 Bool[0, 0, 0, 0], [1, 2, 3, 19], 
 UInt32[0x00000001, 0x00010001, 0x00030000, 0x00000000], 
 3)
        =#
        # clear symbolic hashtable
        # clear matrix
        matrix    = initialize_matrix(ring, C)
        symbol_ht = initialize_secondary_hashtable(hashtable)

        if i > 10_000
            # TODO: log useful info here
            @log level = 100 "Something has probably gone wrong in F4. An exception will follow."
            __error_maximal_number_exceeded(
                "Something has probably gone wrong in F4. Please submit a github issue."
            )
        end
    end

    # remove redundant elements
    filter_redundant!(basis)
    @log "Filtered elements marked redundant"

    @log "Before reduction" basis

    if params.reduced
        @log "Autoreducing the final basis.."
        f4_reducegb_learn!(graph, ring, basis, matrix, hashtable, symbol_ht)
        @log "Autoreduced!"
    end 

    @log "Finalizing computation graph"
    finalize_graph!(graph)

    standardize_basis_learn!(graph, ring, basis, hashtable, hashtable.ord)

    @log "After learn standardization" basis

    #=
┌ Info: After learn standardization
└   basis = Main.Groebner.Basis{UInt64}(
    Vector{Int32}[[4, 5], [6, 3]], 
    Vector{UInt64}[[0x0000000000000001, 0x0000000000000001], 
    [0x0000000000000001, 0x000000007ffffffe]], 
    2, 2, 2, 
    Bool[0, 0], 
    [1, 2], 
    UInt32[0x00000001, 0x00030000], 
    2)
    =#
    # @invariant hashtable_well_formed(:output_f4!, ring, hashtable)
    @invariant basis_well_formed(:output_f4_learn!, ring, basis, hashtable)

    nothing
end

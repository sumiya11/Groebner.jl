
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
    linear_algebra!(graph, ring, matrix, basis, :apply, rng)
    @log "After linear algebra" matrix

    convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht)
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
    nlow = length(nonzeroed_rows)
    nup = matrix_info.nup
    ncols = matrix_info.ncols

    matrix.uprows = Vector{Vector{ColumnIdx}}(undef, nup)
    matrix.lowrows = Vector{Vector{ColumnIdx}}(undef, nlow)
    matrix.low2coef = Vector{Int}(undef, nlow)
    matrix.up2coef = Vector{Int}(undef, nup)

    check_enlarge_hashtable!(symbol_ht, nlow + nup + 2)
    @inbounds for i in 1:nlow
        mult_idx = lowmults[i]
        poly_idx = lowrows[i]

        h = hashtable.hashdata[mult_idx].hash
        etmp = hashtable.monoms[mult_idx]
        rpoly = basis.monoms[poly_idx]

        @log level=1000 "Preserved" poly_idx hashtable.monoms[rpoly]

        @log "Low row" etmp

        # vidx = insert_in_hash_table!(symbol_ht, etmp)

        matrix.lowrows[i] =
            multiplied_poly_to_matrix_row!(symbol_ht, hashtable, h, etmp, rpoly)
        symbol_ht.hashdata[matrix.lowrows[i][1]].idx = 2  # TODO!!

        matrix.low2coef[i] = poly_idx

        @log "low row x2: " matrix.lowrows[i] symbol_ht.monoms[matrix.lowrows[i]]
    end

    @inbounds for i in 1:nup
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

    matrix.ncols = ncols

    matrix.nrows = nlow + nup
    matrix.nlow = nlow
    matrix.nup = nup
    matrix.size = matrix.nrows
end

function f4_apply!(graph, ring, basis::Basis{C}, params) where {M <: Monom, C <: Coeff}
    @invariant basis_well_formed(:input_f4_apply!, ring, basis, graph.hashtable)

    normalize_basis!(ring, basis)

    iters_total = length(graph.matrix_infos)
    iters = 0
    hashtable = graph.hashtable

    matrix = initialize_matrix(ring, C)
    @log "Applying modulo $(ring.ch)"

    update_basis!(basis, hashtable)

    @log "Input basis:" basis
    #=
    ┌ Info: Input basis:
    └   basis = Main.Groebner.Basis{UInt64}(
    Vector{Int32}[[2, 3, 4], [5, 6, 7], [8, 9]], 
    Vector{UInt64}[
    [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], 
    [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], 
    [0x0000000000000001, 0x0000000000000010]], 
    3, 3, 3, Bool[0, 0, 0], [1, 2, 3], 
    UInt32[0x00000001, 0x00000401, 0x00100401], 3)
    =#

    while iters < iters_total
        iters += 1
        @log "F4 Apply iteration $iters"

        symbol_ht = initialize_secondary_hashtable(hashtable)

        symbolic_preprocessing!(graph, iters, basis, matrix, hashtable, symbol_ht)
        @log "After symbolic preprocessing:" matrix
#=
┌ Info: After symbolic preprocessing:
└   matrix = Main.Groebner.MacaulayMatrix{UInt64}(
Vector{Int32}[
[5, 6, 7], 
[6, 7, 8], 
[7, 8, 9], 
[10, 8, 9]], 
Vector{Int32}[[2, 3, 4]], 
Int32[], 
Vector{UInt64}[], 
5, 0, 5, 7, 4, 1, 0, 0, 
[2, 1, 4, 1], 
[4], Int32[], Int32[])
=#

        #=
┌ Info: Apply: after mapping and sorting
└   matrix = Main.Groebner.MacaulayMatrix{UInt64}(
    Vector{Int32}
[[1, 2, 3], 
[2, 3, 6], 
[3, 6, 7], 
[4, 6, 7]], 
Vector{Int32}[[5, 8, 9]], 
Int32[5, 6, 7, 10, 2, 8, 9, 3, 4], 
Vector{UInt64}[], 
5, 0, 5, 9, 4, 1, 5, 4, 
[2, 1, 4, 1], 
[4], 
Int32[], Int32[])
        =#

        reduction_apply!(graph, ring, basis, matrix, hashtable, params.rng, symbol_ht)
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
Vector{Int32}[[2, 3, 4], [5, 6, 7], 
[8, 9], [10, 7, 11], #undef, #undef], 
Vector{UInt64}[[0x0000000000000001, 0x0000000000000001, 0x0000000000000001], [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], [0x0000000000000001, 0x0000000000000010], [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], #undef, #undef], 
6, 4, 4, Bool[0, 0, 0, 0, 0, 0], 
[1, 2, 3, 4, 2, 1438595541200], 
UInt32[0x00000001, 0x00000401, 0x00100401, 0x00000c00, 0x00089805, 0x00000000], 
4)
=#

        matrix = initialize_matrix(ring, C)
    end

    @log "Before reduction" basis
    #=
    ┌ Info: Before reduction
    └   basis = Main.Groebner.Basis{UInt64}(
    Vector{Int32}[[2, 3, 4], [5, 6, 7], [8, 9],
    [10, 7, 11], #undef, #undef], 
    Vector{UInt64}[
    [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], 
    [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], 
    [0x0000000000000001, 0x0000000000000010], 
    [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], 
    #undef, #undef], 
    6, 4, 4, 
    Bool[0, 0, 0, 0, 0, 0], 
    [1, 2, 3, 4, 3, 1437949921248], 
    UInt32[0x00000001, 0x00000401, 0x00100401, 0x00000c00, 0xcedb2f00, 0x0000014e], 
    4)
    =#

    filter_redundant!(basis)

    if params.reduced
        # TODO: drop redundant!!
        @log "Autoreducing the final basis.."
        symbol_ht = initialize_secondary_hashtable(hashtable)
        reducegb_f4!(ring, basis, matrix, hashtable, symbol_ht)
        @log "Autoreduced!"
    end

    standardize_basis!(ring, basis, hashtable, hashtable.ord)

    nothing
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

    @log "Input basis:" basis
    #=
    ┌ Info: Input basis:
    └   basis = Main.Groebner.Basis{UInt64}(Vector{Int32}[
    [2, 3, 4], [5, 6, 7], [8, 9]], 
    Vector{UInt64}[[0x0000000000000001, 0x0000000000000001, 0x0000000000000001], 
    [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], 
    [0x0000000000000001, 0x000000007ffffffe]], 
    3, 3, 3, Bool[0, 0, 0], [1, 2, 3], 
    UInt32[0x00000001, 0x00000401, 0x00100401], 3)
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
Vector{Int32}[
[2, 3, 4], 
[6, 2, 3], 
[3, 4, 8], 
[7, 4, 8]], 
Vector{Int32}[[2, 5], [6, 2, 7]], 
Int32[], Vector{UInt64}[], 
6, 0, 6, 7, 4, 2, 0, 0, 
[1, 2, 4, 1, 4, 1438590237504, 1438504029072, 140732279395472, 9, 10, 11, 12, 13, 14, 15, 16],
[3, 4, 563204, 4], 
Int32[7, 3, 4, 11, 563205, 0, 0, 0, -208472464, 334, -305361008, 334, -264317296, 334, -930281424, 32766], Int32[9, 2, -998074112, 32766])


=#
#=
┌ Info: After symbolic preprocessing:
└   matrix = Main.Groebner.MacaulayMatrix{UInt64}(Vector{Int32}[[5, 6, 7], [6, 7, 8], [7, 8, 9], [10, 8, 9]], Vector{Int32}[[2, 3, 4]], Int32[], Vector{UInt64}[], 5, 0, 5, 7, 4, 1, 0, 0, [2, 1, 4, 1], [4], Int32[], Int32[])
=#

        #=
┌ Info: Learn: after mapping and sorting
└   matrix = Main.Groebner.MacaulayMatrix{UInt64}(
Vector{Int32}[
[1, 2, 3], 
[2, 3, 5], 
[3, 5, 6], 
[4, 5, 6]], 
Vector{Int32}[[2, 7], [1, 2, 4]], 
Int32[6, 2, 3, 7, 4, 8, 5], 
Vector{UInt64}[], 
6, 0, 6, 7, 4, 2, 4, 3, 
[2, 1, 4, 1, 3, 1438590421280, 1438500186192, 140732198948896, 140732279410256, 1437223016464, 140732279394416, 140732279394416, 140732279395296, 140732279395120, 140732279394416, 140732279395120], 
[3, 4, 0, 0], 
Int32[3, 7, 4, 11, -1597329232, 334, 641, 0, 16, 0, -993265568, 32766, -993265568, 32766, -993265568, 32766], 
Int32[9, 2, -1597551784, 334])
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
        Vector{Int32}
        [[1, 3, 4], [2, 4, 5]], 
        Vector{Int32}[[10, 7, 11]], 
        Int32[2, 5, 3, 4, 6], 
        Vector{UInt64}[[0x0000000000000001, 0x0000000000000001, 0x0000000000000001]], 
        1, 1, 1, 5, 2, 1, 2, 3, 
        [1, 1, 563204, 3, 3, 1438140937552, 1438140937584, 1438140937520], 
        [2, 334, 1, 334, -651165439], 
        Int32[3, 4, 0, 0, 563204, 0, 0, 0], 
        Int32[9, 334])

        └   basis = Main.Groebner.Basis{UInt64}(
        Vector{Int32}[[2, 3, 4], [5, 6, 7], 
        [8, 9], [10, 7, 11], #undef, #undef], 
        Vector{UInt64}[[0x0000000000000001, 0x0000000000000001, 0x0000000000000001], [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], [0x0000000000000001, 0x000000007ffffffe], [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], #undef, #undef], 
        6, 3, 4, Bool[0, 0, 0, 0, 0, 0], 
        [1, 2, 3, 4, 4, 1438023461520], 
        UInt32[0x00000001, 0x00000401, 0x00100401, 0x0000014e, 0xa0d23d58, 0x0000014e], 
        3)

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
Vector{Int32}[[2, 3, 4], [5, 6, 7], 
[8, 9], [10, 7, 11], #undef, #undef], 
Vector{UInt64}[[0x0000000000000001, 0x0000000000000001, 0x0000000000000001], [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], [0x0000000000000001, 0x000000007ffffffe], [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], #undef, #undef], 
6, 4, 4, Bool[0, 0, 0, 0, 0, 0], 
[1, 2, 3, 4, 3, 140732277058640], 
UInt32[0x00000001, 0x00000401, 0x00100401, 0x00000c00, 0x00089805, 0x00000000], 
4)
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
    #=
    ┌ Info: Before reduction
    └   basis = Main.Groebner.Basis{UInt64}(
    Vector{Int32}[[2, 3, 4], [5, 6, 7], [8, 9], 
    [10, 7, 11], [13, 9], #undef], 
    Vector{UInt64}[
        [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], 
        [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], 
        [0x0000000000000001, 0x000000007ffffffe], 
        [0x0000000000000001, 0x0000000000000001, 0x0000000000000001], 
        [0x0000000000000001, 0x000000007ffffffe], #undef], 
        6, 5, 5, 
        Bool[0, 0, 0, 0, 0, 0], 
        [1, 2, 3, 4, 5, 0], 
        UInt32[0x00000001, 0x00000401, 0x00100401, 0x00000c00, 0x00700000, 0x00000000], 
        5)
    =#
    if params.reduced
        @log "Autoreducing the final basis.."
        reducegb_f4!(ring, basis, matrix, hashtable, symbol_ht)
        @log "Autoreduced!"
    end

    standardize_basis!(ring, basis, hashtable, hashtable.ord)

    # @invariant hashtable_well_formed(:output_f4!, ring, hashtable)
    # @invariant basis_well_formed(:output_f4!, ring, basis, hashtable)

    nothing
end

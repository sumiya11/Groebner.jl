
#=
    The file contains implementation of the F4 algorithm,
    described originally by Jean-Charles Faugere in
        "A new efficient algorithm for computing Grobner bases"
=#

#------------------------------------------------------------------------------

#

#------------------------------------------------------------------------------

# s-pair, a pair of polynomials
struct SPair
    # first generator as index from the basis array
    poly1::Int
    # second generator -//-
    poly2::Int
    # position of lcm(poly1, poly2) in hashtable
    lcm::Int
    # total degree of lcm
    deg::UInt
end

mutable struct Pairset
    pairs::Vector{SPair}
    # number of filled pairs,
    # Initially zero
    load::Int
end

function initialize_pairset(; initial_size=2^6) # TODO: why 64?
    pairs = Vector{SPair}(undef, initial_size)
    return Pairset(pairs, 0)
end

function Base.isempty(ps::Pairset)
    return ps.load == 0
end

function check_enlarge_pairset!(ps::Pairset, added::Int)
    sz = length(ps.pairs)
    if ps.load + added >= sz
        newsz = max(2 * sz, ps.load + added)
        resize!(ps.pairs, newsz)
    end
end

#------------------------------------------------------------------------------

# TODO: get rid of `Tv` dependency
mutable struct Basis{Tv}
    # vector of polynomials, each polynomial is a vector of monomials,
    # each monomial is represented with it's position in hashtable
    gens::Vector{Vector{Int}}
    # polynomial coefficients
    coeffs::Vector{Vector{Tv}}

    #= Keeping track of sizes   =#
    #=  ndone <= ntotal <= size =#
    # total allocated size,
    # size == length(gens) is always true
    size::Int
    # number of processed polynomials in `gens`
    # Iitially a zero
    ndone::Int
    # total number of polys filled in `gens`
    # (these will be handled during next update! call)
    # Iitially a zero
    ntotal::Int

    #= Keeping track of redundancy =#
    #= invariant =#
    #= length(lead) == length(nonred) == count(isred) == nlead =#
    # if element of the basis
    # is redundant
    isred::Vector{Int8}
    # positions of non-redundant elements in the basis
    nonred::Vector{Int}
    # division masks of leading monomials of
    # non redundant basis elements
    lead::Vector{UInt32}
    # number of filled elements in lead
    nlead::Int
end

function initialize_basis(ring::PolyRing{Tv}, ngens) where {Tv}
    #=
        always true
        length(gens) == length(coeffs) == length(isred) == size
    =#

    sz     = ngens * 2
    ndone  = 0
    ntotal = 0

    gens   = Vector{Vector{Int}}(undef, sz)
    coeffs = Vector{Vector{Tv}}(undef, sz)
    isred  = zeros(Int8, sz)
    nonred = Vector{Int}(undef, sz)
    lead = Vector{UInt32}(undef, sz)

    return Basis{Tv}(gens, coeffs, sz, ndone, ntotal, isred, nonred, lead, 0)
end

function check_enlarge_basis!(basis, added::Int)
    if basis.ndone + added >= basis.size
        basis.size = max(basis.size * 2, basis.ndone + added)
        resize!(basis.gens, basis.size)
        resize!(basis.coeffs, basis.size)
        resize!(basis.isred, basis.size)
        resize!(basis.nonred, basis.size)
        resize!(basis.lead, basis.size)
    end
end

function normalize_basis!(basis)
    cfs = basis.coeffs
    for i in 1:basis.ntotal
        mul = inv(cfs[i][1])
        for j in 2:length(cfs[i])
            cfs[i][j] *= mul
        end
        cfs[i][1] = one(cfs[i][1])
    end
end

#------------------------------------------------------------------------------


function select_normal!(pairset, basis, matrix, ht, symbol_ht; maxpairs=0)

    println("Pairsload: $(pairset.load)")
    println("Critical pairs: $(pairset.pairs[1:pairset.load])")

    sort_pairset_by_degree!(pairset, pairset.load)
    # sort by degree

    ps = pairset.pairs
    min_deg = ps[1].deg

    println("Critical pairs: $(pairset.pairs[1:pairset.load])")
    @info "Min degree is " min_deg

    min_idx = 0
    # beda with boundaries
    while min_idx < pairset.load && ps[min_idx + 1].deg == min_deg
        min_idx += 1
    end

    # number of selected pairs
    npairs = min_idx
    reinitialize_matrix!(matrix, npairs)

    uprows  = matrix.uprows
    lowrows = matrix.lowrows

    # polynomials from pairs in order (p11, p12)(p21, p21)
    # (future rows of the matrix)
    gens   = Vector{Int}(undef, 2*npairs)

    # buffer !
    etmp = ht.exponents[1]

    @info "Other things " min_idx npairs

    i = 1
    while i <= npairs
        matrix.ncols += 1
        load = 1
        lcm = ps[i].lcm
        j = i

        @info "handling pair $i with lcm $lcm"
        @info "btw, its exponent is" ht.exponents[lcm]
        # we collect all generators with same lcm into gens

        while j <= npairs && ps[j].lcm == lcm
            gens[load] = ps[j].poly1
            load += 1
            gens[load] = ps[j].poly2
            load += 1
            j += 1
        end
        load -= 1

        @info "collected $(load) polys with same lcm"
        println("gens now:", gens)


        # sort by sparsity ?
        # upd: yes!
        # upd: no!
        # sort by number in the basis (sort by itself)
        sort_generators_by_something!(gens, load, basis)

        println("sorted gens:", gens)

        # now we collect reducers, and reduced

        # first generator index in groebner basis
        prev = gens[1]
        # first generator in hash table
        poly = basis.gens[prev]
        # first generator lead monomial index in hash data
        vidx = poly[1]

        # first generator exponent
        eidx = ht.exponents[vidx]
        # exponent of lcm corresponding to first generator
        elcm = ht.exponents[lcm]

        @info "? divides" elcm eidx

        for u in 1:ht.explen
            etmp[u] = elcm[u] - eidx[u]
        end
        # now etmp contents complement to eidx in elcm
        @info "q = " etmp

        # hash of complement
        htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

        # add row as a reducer
        matrix.nup += 1
        uprows[matrix.nup] = multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly)
        # map upper row to index in basis
        matrix.up2coef[matrix.nup] = prev
        @info "inserted reducer row " matrix.nup uprows[matrix.nup]

        # TODO: nreduce increment

        # mark lcm column as reducer in symbolic hashtable
        symbol_ht.hashdata[uprows[matrix.nup][1]].idx = 2
        matrix.nrows += 1

        # over all polys with same lcm
        for k in 1:load
            # duplicate generator
            if gens[k] == prev
                continue
            end

            # if the table was reallocated
            elcm = ht.exponents[lcm]

            # index in gb
            prev = gens[k]
            # poly of indices of monoms in hash table
            poly = basis.gens[prev]
            vidx = poly[1]
            # leading monom idx
            eidx = ht.exponents[vidx]
            for u in 1:ht.explen
                etmp[u] = elcm[u] - eidx[u]
            end

            htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

            # add row to be reduced
            matrix.nlow += 1
            lowrows[matrix.nlow] = multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly)
            # map lower row to index in basis
            matrix.low2coef[matrix.nlow] = prev
            @info "inserted reduced row " matrix.nlow lowrows[matrix.nlow]

            symbol_ht.hashdata[lowrows[matrix.nlow][1]].idx = 2

            matrix.nrows += 1

        end

        i = j
    end

    resize!(matrix.lowrows, matrix.nrows - matrix.ncols)

    # remove selected parirs from pairset
    for i in 1:pairset.load - npairs
        ps[i] = ps[i + npairs]
    end
    pairset.load -= npairs

end

#------------------------------------------------------------------------------

# Prepares the given set of polynomials L to be reduced by the ideal G,
# starts the reduction routine,
# and returns reduced polynomials
function reduction!(basis, matrix, ht, symbol_ht; linalg=:sparse)

     convert_hashes_to_columns!(matrix, symbol_ht)

     sort_matrix_rows_decreasing!(matrix) # for pivots,  A part
     sort_matrix_rows_increasing!(matrix) # for reduced, B part

     printstyled("### PREPARED MATRIX ###\n", color=:red)
     # matrix.ncols beda

     dump(matrix, maxdepth=2)

     exact_sparse_linear_algebra!(matrix, basis)


     printstyled("### AFTER RREF ###\n", color=:red)
     dump(matrix, maxdepth=4)
     
     convert_matrix_rows_to_basis_elements!(matrix, basis, ht, symbol_ht)
end

#------------------------------------------------------------------------------

function find_multiplied_reducer!(basis, matrix, ht, symbol_ht, vidx)

    e = symbol_ht.exponents[vidx]
    etmp = ht.exponents[1]
    divmask = symbol_ht.hashdata[vidx].divmask

    blen = basis.ndone     #
    leaddiv = basis.lead  #

    # searching for a poly from state.basis whose leading monom
    # divides the given exponent e
    i = 1
    @label Letsgo
    while i <= blen && (leaddiv[i] & divmask) != 0
        i += 1
    end

    # here found polynomial from basis with leading monom
    # dividing symbol_ht.exponents[vidx]
    if i <= blen
        # reducers index and exponent in hash table
        rpoly = basis.gens[i]
        rexp  = ht.exponents[rpoly[1]]

        for j in 1:ht.explen
            # if it actually does not divide and divmask lies
            if e[j] < rexp[j]
                i += 1
                @goto Letsgo
            end
            etmp[j] = e[j] - rexp[j]
        end
        # now etmp = e // rexp in terms of monomias,
        # hash is linear
        h = symbol_ht.hashdata[vidx].hash - ht.hashdata[rpoly[1]].hash

        nreducer = matrix.nup
        matrix.uprows[matrix.nup + 1] = multiplied_poly_to_matrix_row!(symbol_ht, ht, h, copy(etmp), rpoly)
        matrix.up2coef[matrix.nup + 1] = rpoly

        # upsize matrix
        symbol_ht.hashdata[vidx].idx = 2
        matrix.nup += 1
        nreducer += 1
    end

end

# Given the set of polynomials L and the basis G,
# extends L to contain possible polys for reduction by G,
# and returns it
function symbolic_preprocessing!(
                basis, matrix,
                ht, symbol_ht;
                linalg=:sparse)

    # check matrix sizes (we want to omit this I guess)

    symbol_load = symbol_ht.load

    nrr    = matrix.ncols
    onrr   = matrix.ncols

    @info "Entering symbolic symbolic_preprocessing! with " symbol_load nrr onrr

    #=
        TODO: matrix_enlarge!
    =#
    # wait a minute..
    while matrix.size <= nrr + symbol_load
        matrix.size *= 2
        resize!(matrix.uprows, matrix.size)
    end

    @info "Now matrix is of size $(matrix.size)"

    # for each lcm present in symbolic_ht set on select stage
    # !!! remember that first index is for buffer
    # TODO
    i = 2
    #= First round, we add multiplied polynomials which divide  =#
    #= a monomial exponent from selected spairs                 =#
    while i <= symbol_load
        # not a reducer
        if symbol_ht.hashdata[i].idx == 0
            @info "find_multiplied_reducer for $i"
            symbol_ht.hashdata[i].idx = 1
            matrix.ncols += 1
            find_multiplied_reducer!(basis, matrix, ht, symbol_ht, i)
        end
        i += 1
    end

    @info "round 2"

    #= Second round, we add multiplied polynomials which divide  =#
    #= lcm added on previous for loop                             =#
    while i <= symbol_ht.load
        if matrix.size == nrr
            matrix.size *= 2
            resize!(matrix.uprows, matrix.size)
        end
        @info "find_multiplied_reducer for $i"

        symbol_ht.hashdata[i].idx = 1
        matrix.ncols += 1
        find_multiplied_reducer!(basis, matrix, ht, symbol_ht, i)
        i += 1
    end

    # shrink matrix sizes, set constants
    resize!(matrix.uprows, nrr)
    matrix.nrows += nrr - onrr
    matrix.nup   = nrr
    matrix.nlow  = matrix.nrows - nrr
    matrix.size  = matrix.nrows

end



function update_pairset!(
            pairset,
            basis,
            ht,
            update_ht,
            idx)

    pl = pairset.load
    bl = basis.ntotal
    nl = pl + bl
    ps = pairset.pairs

    new_lead  = basis.gens[idx][1]

    @info "entering update_pairset" pl bl nl ps
    @info "new poly is " new_lead ht.exponents[new_lead]

    # initialize new critical lcms
    plcm = Vector{Int}(undef, basis.ntotal + 1)

    # TODO
    reinitialize_hash_table!(update_ht, basis.ntotal)

    # for each combination (new_Lead, basis.gens[i][1])
    # generate a pair
    for i in 1:basis.ntotal
        plcm[i] = get_lcm(basis.gens[i][1], new_lead, ht, update_ht)
        @debug "new plcm" plcm[i]
        deg = update_ht.hashdata[plcm[i]].deg
        newidx = pl + i
        if basis.isred[i] == 0
            ps[newidx] = SPair(i, idx, plcm[i], deg)
        else
            # lcm == 0 will mark redundancy of spair
            ps[newidx] = SPair(i, idx, 0, deg)
        end
    end

    @debug "INSIDE 1" ht update_ht

    @info "generated pairset" pairset

    # traverse existing pairs
    for i in 1:pl
        j = ps[i].poly1
        l = ps[i].poly2
        m = max(ps[pl + l].deg, ps[pl + j].deg)

        @info "traverse existing" j l m
        # if an existing pair is divisible by lead of new poly
        # and has a greater degree than newly generated one
        if is_monom_divisible(ps[i].lcm, new_lead, ht) && ps[i].deg > m
            # mark lcm as 0
            ps[i] = SPair(ps[i].poly1, ps[i].poly2, 0, ps[i].deg)
        end
    end
    # TODO: this can be done faster is we
    # create only non-redundant pairs at first

    # traverse new pairs to check for redundancy
    j = 1
    for i in 1:bl
        if basis.isred[i] == 0
            ps[pl + j] = ps[pl + i]
            j += 1
        end
    end
    @info "after sifting through new pairs" pairset j

    sort_pairset_by_degree!(pairset, j - 1)
    @info "after sorting" pairset

    for i in 1:j - 1
        plcm[i] = ps[pl + i].lcm
    end
    plcm[j] = 0
    pc = j

    @info "plcm" plcm pc
    pc -= 1

    # mark redundancy of some pairs from plcm
    for j in 1:pc
        # if is not redundant
        if plcm[j] > 0
            check_monomial_division_in_update(plcm, j, pc, plcm[j], update_ht)
        end
    end

    # remove useless pairs from pairset
    # by moving them to the end
    j = 1
    for i in 1:pl
        ps[i].lcm == 0 && continue
        ps[j] = ps[i]
        j += 1
    end

    @info "after removing redundant" ps

    # assure that basis hashtable can store new lcms
    if ht.size - ht.load <= pc
        enlarge_hash_table!(ht)
    end

    # add new lcms to the basis hashtable
    insert_plcms_in_basis_hash_table!(pairset, pl, ht, update_ht, basis, plcm, j, pc)

    @info "after insert plcms" ht

    # mark redundant elements in masis
    nonred = basis.nonred
    lml = basis.nlead
    for i in 1:lml
        if basis.isred[nonred[i]] == 0
            if is_monom_divisible(basis.gens[nonred[i]][1], new_lead, ht)
                basis.red[nonred[i]] = 1
            end
        end
    end

end

function update_basis!(
            basis,
            ht::CustomMonomialHashtable,
            update_ht::CustomMonomialHashtable)

    # here we could check overall redundancy and update basis.lead

    k = 1
    lead   = basis.lead
    nonred = basis.nonred

    @info "basis" basis.ndone basis.ntotal
    @info "updating basis & lead" lead nonred basis.isred basis.nlead
    for i in 1:basis.nlead
        if basis.isred[nonred[i]] == 0
            basis.lead[k]   = lead[i]
            basis.nonred[k] = nonred[i]
            k += 1
        end
    end
    basis.nlead = k - 1

    @info "updated #1" lead nonred basis.isred basis.nlead

    for i in basis.ndone+1:basis.ntotal
        if basis.isred[i] == 0
            lead[k]   = ht.hashdata[basis.gens[i][1]].divmask
            nonred[k] = i
            k += 1
        end
    end

    basis.nlead = k - 1
    basis.ndone = basis.ntotal

    @info "updated #2" lead nonred basis.isred basis.nlead
end

# checks if element of basis at position idx is redundant
function is_redundant!(
            pairset, basis, ht, update_ht, idx)

    pairs = pairset.pairs

    # lead of new polynomial
    lead_new = basis.gens[idx][1]
    # degree of lead
    lead_deg = ht.hashdata[lead_new].deg

    start = basis.ndone + 1
    for i in 1:idx - 1
        if basis.isred[i] == 1
            continue
        end

        lead_i = basis.gens[i][1]

        @info "?? divisibility" lead_new lead_i
        @info "" ht.exponents[lead_new] ht.exponents[lead_i]

        if is_monom_divisible(lead_new, lead_i, ht)
            println("Delit !")
            lcm_new = get_lcm(lead_i, lead_new, ht, ht)
            @info "lcm is " lcm_new ht.exponents[lcm_new]

            psidx = pairset.load + 1
            pairs[psidx] = SPair(i, idx, lcm_new, ht.hashdata[lcm_new].deg)

            basis.isred[idx] = 1
            pairset.load += 1

            return true
        end
    end

    return false
end

# I am doing this

function update!(
        pairset::Pairset,
        basis::Basis{Tv},
        ht::CustomMonomialHashtable,
        update_ht::CustomMonomialHashtable) where {Tv}

    # reinitialize_hash_table!(update_ht, length(basis.gens))

    # number of added elements
    npivs = basis.ntotal - basis.ndone

    # number of potential critical pairs to add
    npairs = basis.ndone * npivs
    for i in 1:npivs
        npairs += i
    end

    @info "" basis.ndone basis.ntotal
    @info "potentially adding $npairs pairs from $npivs pivots"

    # make sure pairset and update hashtable have enough
    # space to store new pairs
    check_enlarge_pairset!(pairset, npairs)

    # TODO: better to reinitialize here, or in is_redundant! ?
    reinitialize_hash_table!(update_ht, basis.ntotal)

    @debug "status ONE" pairset ht update_ht

    # for each new element in basis
    for i in basis.ndone+1:basis.ntotal
        @info "adding $i th poly in update"
        # check redundancy of new poly
        if is_redundant!(pairset, basis, ht, update_ht, i)
            @info "Redundant!"
            continue
        end
        @info "Not redundant!"

        @debug "status TWO $i" pairset ht update_ht

        update_pairset!(pairset, basis, ht, update_ht, i)

        @debug "status THREE $i" pairset ht update_ht
    end

    update_basis!(basis, ht, update_ht)
end

#------------------------------------------------------------------------------

# given input exponent and coefficient vectors hashes exponents into `ht`
# and then constructs hashed polynomial vectors for `basis`
function fill_data!(basis, ht, exponents, coeffs)
    etmp = ht.exponents[1]
    ngens = length(exponents)

    for i in 1:ngens
        while length(exponents[i]) >= ht.size - ht.load
            enlarge_hash_table!(ht)
        end

        nterms = length(coeffs[i])
        basis.coeffs[i] = coeffs[i]
        basis.gens[i] = Vector{Int}(undef, nterms)
        poly = basis.gens[i]
        for j in 1:nterms
            poly[j] = insert_in_hash_table!(ht, exponents[i][j])
        end

        # sort terms (not needed),
        # beautify coefficients (not needed),
        # and all
    end

    basis.ntotal = ngens
end

function initialize_structures(
            ring::PolyRing,
            exponents::Vector{Vector{Vector{UInt16}}},
            coeffs::Vector{Vector{Tv}},
            tablesize) where {Tv}


    # basis for storing basis elements,
    # pairset for storing critical pairs of basis elements to assess,
    # hashtable for hashing monomials occuring in the basis
    basis    = initialize_basis(ring, length(exponents))
    basis_ht = initialize_basis_hash_table(ring)

    # filling the basis and hashtable with the given inputs
    fill_data!(basis, basis_ht, exponents, coeffs)

    # every monomial in hashtable is associated with its divmask
    # to perform divisions faster. Filling those
    fill_divmask!(basis_ht)

    # sort input, smaller leading terms first
    sort_gens_by_lead_increasing!(basis, basis_ht)

    @info "sorted basis gens" basis.gens
    @info "leads:" [basis_ht.exponents[e] for e in first.(basis.gens[1:basis.ntotal])]

    # divide each polynomial by leading coefficient
    normalize_basis!(basis)

    basis, basis_ht
end

#------------------------------------------------------------------------------

function dump_all(msg, pairset, basis, matrix, ht, update_ht, symbol_ht)
    printstyled("### $msg ###\n", color=:yellow)
    printstyled("internal repr\n", color=:red)

    printstyled(". pairset filled $(pairset.load) / $(length(pairset.pairs))\n", color=:red)
    println(pairset.pairs[1:pairset.load+1])

    printstyled(". basis filled $(basis.ntotal) / $(length(basis.gens)), where $(basis.ndone) are done\n", color=:red)
    dump(basis, maxdepth=3)

    printstyled(". matrix size $(matrix.size), npivs $(matrix.npivots), nrows $(matrix.nrows), ncols $(matrix.ncols)\n", color=:red)
    printstyled("nup $(matrix.nup), nlow $(matrix.nlow), nleft $(matrix.nleft), nright $(matrix.nright)\n", color=:red)
    dump(matrix, maxdepth=2)

    printstyled(". basis hashtable filled $(ht.load) / $(ht.size)\n", color=:red)
    dump(ht, maxdepth=3)

    printstyled(". update hashtable filled $(update_ht.load) / $(update_ht.size)\n", color=:red)
    dump(update_ht, maxdepth=2)

    printstyled(". symbol hashtable filled $(symbol_ht.load) / $(symbol_ht.size)\n", color=:red)
    dump(symbol_ht, maxdepth=2)

    printstyled("nice repr\n", color=:green)



    printstyled("##################\n", color=:yellow)
end

#------------------------------------------------------------------------------

function f4(polys::Vector{MPoly{Tv}};
            select=select_normal!,
            reduced=true,
            linalg=:sparse,
            maxpairs=0,
            tablesize=2^16) where {Tv}

    ring, exps, coeffs = convert_to_internal(polys)
    gb, ht = f4(ring, exps, coeffs;
        select=select, reduced=reduced,
        linalg=linalg, maxpairs=maxpairs, tablesize=tablesize )

    export_basis(parent(first(polys)), gb, ht)
end

#
"""
    Main function to calculate the Groebner basis of the given polynomial ideal.
    Specialized to work only over finite fields.

    Parameters
        . F         - an array of polynomials over finite field
        . select    - a strategy for polynomial selection on each iteration
        . reduced   - reduce the basis so that the result is unique
        . linalg    - linear algebra backend to use. Possible options
                      are :dense , :sparse , and :sparserand. Default is :sparse
        . maxpairs  - maximal number of pairs selected for one matrix; default is
                      0, i.e. no restriction. If matrices get too big or consume
                      too much memory this is a good parameter to play with.
        . tablesize -
"""
function f4(ring::PolyRing,
            exponents::Vector{Vector{Vector{UInt16}}},
            coeffs::Vector{Vector{Tv}};
            select=selectnormal,
            reduced=true,
            linalg=:sparse,
            maxpairs=0,
            tablesize=2^16) where {Tv}

    printstyled("### INPUT ###\n", color=:red)
    println("coeffs = $coeffs")
    println("exps = $exponents")

    # initialize basis and hashtable with input polynomials,
    # fields are not meaningful yet and will be set during update! step
    basis, ht = initialize_structures(ring, exponents, coeffs, tablesize)

    # matrix storing coefficients in rows
    # wrt columns representing the current monomial basis
    matrix = initialize_matrix(ring)

    # initialize hash tables for update and symbolic preprocessing steps
    update_ht  = initialize_secondary_hash_table(ht)
    symbol_ht  = initialize_secondary_hash_table(ht)

    # a set to store critical pairs of polynomials to be reduced
    pairset = initialize_pairset()

    dump_all("AFTER INIT", pairset, basis, matrix, ht, update_ht, symbol_ht)


    # makes basis fields valid,
    # does not copy,
    # checks for redundancy of new elems
    update!(pairset, basis, ht, update_ht)
    dump_all("INITIAL UPDATE", pairset, basis, matrix, ht, update_ht, symbol_ht)

    d = 0
    # while there are pairs to be reduced
    while !isempty(pairset)
        d += 1
        @warn "F4 ITER $d"

        # selects pairs for reduction from pairset following normal strategy
        # (minimal lcm degrees are selected),
        # and puts these into the matrix rows
        select_normal!(pairset, basis, matrix, ht, symbol_ht, maxpairs=maxpairs)
        dump_all("AFTER SELECT", pairset, basis, matrix, ht, update_ht, symbol_ht)

        symbolic_preprocessing!(basis, matrix, ht, symbol_ht)
        dump_all("AFTER SYMBOLIC", pairset, basis, matrix, ht, update_ht, symbol_ht)

        # reduces polys and obtains new potential basis elements
        reduction!(basis, matrix, ht, symbol_ht, linalg=linalg)
        dump_all("AFTER REDUCTION", pairset, basis, matrix, ht, update_ht, symbol_ht)

        # update the current basis with polynomials produced from reduction,
        # does not copy,
        # checks for redundancy
        update!(pairset, basis, ht, update_ht)
        dump_all("AFTER UPDATE", pairset, basis, matrix, ht, update_ht, symbol_ht)

        if d > 1000
            @error "Something is probably wrong in f4.."
            break
        end
    end

    if reduced
        # TODO
        # reduce_basis!(basis, matrix)
    end

    basis, ht
end

#------------------------------------------------------------------------------

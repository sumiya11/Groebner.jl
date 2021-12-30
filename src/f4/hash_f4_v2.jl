
#=
    The file contains implementation of the F4 algorithm,
    described originally by Jean-Charles Faugere in
        A new efficient algorithm for computing Grobner bases
=#

#------------------------------------------------------------------------------


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
    load::Int
end

function initialize_pairset()
    pairs = Vector{SPair}(undef, 64) # TODO: why 64?
    return Pairset(pairs, 0)
end

function Base.isempty(ps::Pairset)
    return isempty(ps.pairs)
end

function check_enlarge_pairset!(ps::Pairset, added::Int)
    sz = length(ps.pairs)
    if ps.load + added >= sz
        newsz = max(2 * sz, ps.load + added)
        resize!(ps.pairs, newsz)
    end
end

#------------------------------------------------------------------------------

mutable struct Basis{Tv}
    # vector of polynomials, each polynomial is a vector of monomials,
    # each monomial is represented with it's position in hashtable
    gens::Vector{Vector{Int}}
    # leading monomials of gens as short divmasks
    lead::Vector{UInt32}
    # polynomial coefficients
    coeffs::Vector{Vector{Tv}}

    # if element of the basis
    # is redundant
    redundant::Vector{Int8}

    # allocated size
    size::Int

    # number of processed polynomials in `gens`
    # Iitially a zero
    ndone::Int

    # total number of polys filled in `gens`
    # (these will be handled during next update! call)
    # Iitially a zero
    ntotal::Int

    #=
        ndone <= ntotal <= size
    =#
end

function initialize_basis(ring::PolyRing{Tv}, ngens) where {Tv}
    sz     = ngens * 2
    ndone  = 0
    ntotal = 0
    gens   = Vector{Vector{Int}}(undef, ngens)
    lead   = Vector{UInt32}(undef, ngens)
    coeffs = Vector{Vector{Tv}}(undef, ngens)
    red    = zeros(Int8, ngens)
    return Basis{Tv}(gens, lead, coeffs, red, 0, 0)
end

function normalize_basis!(basis)
    cfs = basis.coeffs
    for i in 1:length(cfs)
        cf = cfs[i][1]
        for j in 1:length(cfs[i])
            cfs[i][j] //= cf
        end
    end
end

#------------------------------------------------------------------------------

# Normal selection strategy
# Given an array of pairs, selects all of the lowest degree
# where the degree of (f, g) is equal to the degree of lcm(lm(f), lm(g))
function select_normal!(pairset, basis, matrix, ht, symbol_ht; maxpairs=0)

    println("Pairsload: $(pairset.load)")
    println("Critical pairs: $(pairset.pairs[1:pairset.load])")

    sort_pairset_by_degree!(pairset)
    # sort by degree

    pairs = pairset.pairs
    min_deg = pairs[1].deg

    println("Critical pairs: $(pairset.pairs[1:pairset.load])")
    @info "Min degree is " min_deg

    min_idx = 0
    # beda with boundaries
    while min_idx < pairset.load && pairs[min_idx + 1].deg == min_deg
        min_idx += 1
    end

    # number of selected pairs
    npairs = min_idx
    reinitialize_matrix!(matrix, npairs)

    reducer = matrix.reducer
    reduced = matrix.reduced

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
        lcm = pairs[i].lcm
        j = i

        @info "handling pair $i with lcm $lcm"
        @info "btw, its exponent is" ht.exponents[lcm]
        # we collect all generators with same lcm into gens

        while j <= npairs && pairs[j].lcm == lcm
            gens[load] = pairs[j].poly1
            load += 1
            gens[load] = pairs[j].poly2
            load += 1
            j += 1
        end
        load -= 1

        @info "collected $(load) polys with same lcm"
        println("gens now:", gens)


        # sort by sparsity ?
        # upd: yes!
        # upd: no!
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
        reducer[matrix.nup] = multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly)
        @info "inserted reducer row " matrix.nup reducer[matrix.nup]

        # TODO: nreduce increment

        # mark lcm column as reducer in symbolic hashtable
        symbol_ht.hashdata[reducer[matrix.nup][1]].idx = 2
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
            # leading monom idx
            eidx = ht.exponents[poly[1]]
            for u in 1:ht.explen
                etmp[u] = elcm[u] - eidx[u]
            end

            htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

            # add row to be reduced
            matrix.nlow += 1
            reduced[matrix.nlow] = multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly)

            @info "inserted reduced row " matrix.nlow reduced[matrix.nlow]

            symbol_ht.hashdata[reduced[matrix.nlow][1]].idx = 2

            matrix.nrows += 1

        end

        i = j
    end

    resize!(matrix.reduced, matrix.nrows - matrix.ncols)

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

     convert_matrix_rows_to_basis_elems!(matrix, basis, symbol_ht)
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
        matrix.reducer[matrix.nup + 1] = multiplied_poly_to_matrix_row!(symbol_ht, ht, h, etmp, rpoly)

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
    while matrix.size <= nrr + symbol_load
        matrix.size *= 2
        resize!(matrix.reducer, matrix.size)
    end

    @info "Now matrix is of size $(matrix.size)"

    # for each lcm present in symbolic_ht set on select stage
    # !!! remember that first index is for buffer
    # TODO
    i = 2
    #= First round, we add multiplied polynomials which divide  =#
    #= a monomial exponent from selected spairs                 =#
    while i < symbol_load
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
    while i < symbol_ht.load
        if matrix.size == nrr
            matrix.size *= 2
            resize!(matrix.reducer, matrix.size)
        end
        @info "find_multiplied_reducer for $i"

        symbol_ht.hashdata[i].idx = 1
        matrix.ncols += 1
        find_multiplied_reducer!(basis, matrix, ht, symbol_ht, i)
        i += 1
    end

    # shrink matrix sizes, set constants
    resize!(matrix.reducer, nrr)
    matrix.nrows += nrr - onrr
    matrix.nup   = nrr
    matrix.nlow  = matrix.nrows - nrr
    matrix.size  = matrix.nrows

end

#------------------------------------------------------------------------------

# Hereinafter a set of heuristics is defined to be used in update! (see below)
# They assess which polynomials are worthy of adding to the current pairset

function update_heuristic_1!(
                basis_candidates,
                state,
                newpoly,
                basis_ht, update_ht)

    lead_new = leading_monomial(newpoly)

    for i in 1:state.blength
        lead_i = state.lead[i]
        if is_monom_gcd_constant(lead_i, lead_new, basis_ht)
            basis_candidates[i] = 1
            continue
        end

        divflag = true
        for j in i+1:state.blength
            lead_j = state.lead[j]
            if is_monom_divisible(
                        monom_lcm(lead_i, lead_new, basis_ht, update_ht),
                        monom_lcm(lead_j, lead_new, basis_ht, update_ht),
                        basis_ht)
                divflag = false
                break
            end
        end
        if divflag
            basis_candidates[i] = 1
            continue
        end

        divflag = true
        for j in 1:i-1
            basis_candidates[j] == 0 && continue
            lead_j = state.lead[j]
            if is_monom_divisible(
                        monom_lcm(lead_i, lead_new, basis_ht, update_ht),
                        monom_lcm(lead_j, lead_new, basis_ht, update_ht),
                        basis_ht)
                divflag = false
                break
            end
        end
        if divflag
            basis_candidates[i] = 1
        end
    end
end

function update_heuristic_2!(
                basis_candidates,
                state,
                newpoly,
                basis_ht, update_ht)

    lead_new = leading_monomial(newpoly)

    for i in 1:state.blength
        basis_candidates[i] == 0 && continue
        lead_i = state.lead[i]
        if is_monom_gcd_constant(lead_i, lead_new, basis_ht)
            basis_candidates[i] = 0
        end
    end
end

function update_heuristic_3!(
                basis_candidates,
                state,
                newpoly,
                basis_ht, update_ht)

    lead_new = leading_monomial(newpoly)

    for i in 1:state.blength
        basis_candidates[i] == 0 && continue
        lead_i = state.lead[i]

        lcm_inew = monom_lcm(lead_i, lead_new, basis_ht, update_ht)
        if !is_monom_divisible(lcm_inew, lead_new, basis_ht)
            basis_candidates[i] = 1
            continue
        end


    end
end

function update_heuristic_4!(G, h, ht)
    filter!(g -> !is_monom_divisible(leading_monomial(g), leading_monomial(h), ht), G)
end

function update_pairset!(
            pairset,
            basis,
            ht,
            update_ht,
            idx)

    new_lcm  = basis.gens[idx][1]
    critical = Vector{UInt}(undef, basis.ndone)

    # initialize new critical pairs
    pload = pairset.load

    for i in 1:basis.ndone
        critical[i] = get_lcm(basis.gens[i][1], new_lcm, ht, update_ht)
        deg = update_ht.hashdata[critical[i]].deg
        if basis.redundant[i] == 0
            pairset.pairs[pload] = SPair(i, idx, critical[i], deg)
            pload ++ 1
        end
    end

    basis.ndone  += 1
    # pairset.load += length(critical)

    # TODO
    # update_heuristic_1!(basis_candidates, state, newpoly, basis_ht, update_ht)
    # update_heuristic_2!(basis_candidates, state, newpoly, basis_ht, update_ht)
    # update_heuristic_3!(basis_candidates, state, newpoly, basis_ht, update_ht)

    if ht.size - ht.load <= length(critical)
        enlarge_hash_table!(ht)
    end

    insert_plcms_in_basis_hash_table!(
        critical, ht, update_ht, basis,
        1, length(critical))

    # mark redundant (no redundant for now)
end

function update_basis!(
            basis,
            basis_ht::CustomMonomialHashtable,
            update_ht::CustomMonomialHashtable)

    # here we could check overall redundancy and update basis.lead

end

# checks if element of basis at position idx is redundant
function is_redundant!(
            pairset, basis, ht, update_ht, idx)

    lead_new = basis.gens[idx][1]
    pairs    = pairset.pairs
    for i in 1:basis.ndone
        if basis.redundant[i] == 1
            continue
        end
        lead_i = basis.gens[i][1]
        @info "?? divisibility" lead_new lead_i
        @info "" ht.exponents[lead_new] ht.exponents[lead_i]
        if is_monom_divisible(lead_new, lead_i, ht)
            println("Delit !")
            lcm_new = get_lcm(lead_i, lead_new, ht, ht)
            @info "lcm is " lcm_new ht.exponents[lcm_new]

            pairs[pairset.load+1] = SPair(i, idx, lcm_new, ht.hashdata[lcm_new].deg)
            basis.redundant[idx] = 1
            basis.ndone  += 1
            pairset.load += 1
            return true
        end
    end

    return false
end

# I am doing this

# "Adds" h to the set of generators G and set of pairs P, while applying some
# heuristics to reduce the number of pairs on fly
# Returns new generator and pair sets G, P
function update!(
        pairset, basis, ht, update_ht) where {Tv}


    # reinitialize_hash_table!(update_ht, length(basis.gens))

    # number of added elements
    npivots = basis.ntotal - basis.ndone
    # number of potential critical pairs to add
    npairs = basis.ndone * npivots
    for i in 1:npivots
        npairs += i
    end

    @info "" basis.ndone basis.ntotal
    @info "potentially adding $npairs pairs from $npivots pivots"

    reinitialize_hash_table!(update_ht, basis.ntotal)

    # for each new element in basis
    for i in basis.ndone+1:basis.ntotal
        @info "adding $i th poly in update"
        # check redundancy of new poly
        if is_redundant!(pairset, basis, ht, update_ht, i)
            continue
        end

        @info "Not redundant!"

        update_pairset!(pairset, basis, ht, update_ht, i)
        # TODO: add pairs for basis_candidates

    end

    update_basis!(basis, ht, update_ht)
end

#------------------------------------------------------------------------------

function fill_data!(basis, ht, exponents, coeffs)
    etmp = ht.exponents[1]
    ngens = length(exponents)

    for i in 1:ngens
        while length(exponents[i]) >= ht.size - ht.load
            enlarge_hash_table!(ht)
            etmp = ht.exponents[1]
        end

        nterms = length(coeffs[i])
        basis.coeffs[i] = coeffs[i]
        basis.gens[i]   = Vector{Int}(undef, nterms)
        poly = basis.gens[i]
        for j in 1:nterms
            poly[j] = insert_in_hash_table!(ht, exponents[i][j])
        end

        # sort terms,
        # beautify coefficients,
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
    pairset  = initialize_pairset()
    basis_ht = initialize_basis_hash_table(ring)

    # filling the basis and hashtable with the given inputs
    fill_data!(basis, basis_ht, exponents, coeffs)

    # every monomial in hashtable is associated with its divmask
    # to perform divisions faster. Filling those
    fill_divmask!(basis_ht)

    # divide each polynomial by leading coefficient
    normalize_basis!(basis)

    basis, pairset, basis_ht
end

#------------------------------------------------------------------------------

function dump_all(msg, pairset, basis, matrix, ht, update_ht, symbol_ht)
    printstyled("### $msg ###\n", color=:yellow)
    printstyled("internal repr\n", color=:red)

    printstyled(". pairset filled $(pairset.load) / $(length(pairset.pairs))\n", color=:red)
    println(pairset.pairs[1:pairset.load+1])

    printstyled(". basis filled $(basis.ntotal) / $(length(basis.gens)), where $(basis.ndone) are done\n", color=:red)
    dump(basis, maxdepth=2)

    printstyled(". matrix size $(matrix.size), npivs $(matrix.npivots), nrows $(matrix.nrows), ncols $(matrix.ncols)\n", color=:red)
    printstyled("nup $(matrix.nup), nlow $(matrix.nlow), nleft $(matrix.nleft), nright $(matrix.nright)\n", color=:red)
    dump(matrix, maxdepth=2)

    printstyled(". basis hashtable filled $(ht.load) / $(ht.size)\n", color=:red)
    dump(ht, maxdepth=2)

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
    f4(ring, exps, coeffs;
        select=select, reduced=reduced,
        linalg=linalg, maxpairs=maxpairs, tablesize=tablesize )
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
    basis, pairset, ht = initialize_structures(ring, exponents, coeffs, tablesize)

    # matrix storing coefficients in rows
    # wrt columns representing the current monomial basis
    matrix = initialize_matrix()

    # initialize hash tables for update and symbolic preprocessing steps
    update_ht  = initialize_secondary_hash_table(ht)
    symbol_ht  = initialize_secondary_hash_table(ht)

    # a set to store pairs of polynomials to be reduced
    pairset = initialize_pairset()

    dump_all("AFTER INIT", pairset, basis, matrix, ht, update_ht, symbol_ht)


    # makes basis fields valid, creates
    # does not copy,
    # checks for redundancy
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
        reduce_basis!(basis, matrix)
    end

    basis
end

#------------------------------------------------------------------------------

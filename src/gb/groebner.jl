
function remove_zeros_from_input!(ring::PolyRing, 
            exps::Vector{Vector{M}}, 
            coeffs::Vector{Vector{T}}) where {M, T}
    @assert length(exps) == length(coeffs)
    filter!(!iszero_coeffvector, coeffs)
    filter!(!iszero_monomvector, exps)
    @assert length(exps) == length(coeffs)
    iszerobasis = isempty(exps)
    if iszerobasis
        push!(exps, Vector{M}())
        push!(coeffs, Vector{T}())
    end
    iszerobasis
end

#------------------------------------------------------------------------------

function groebner(polynomials, representation, reduced, ordering, certify, forsolve, linalg, rng)
    #= extract ring information, exponents and coefficients
       from input polynomials =#
    # Copies input, so that polynomials would not be changed itself.
    ring, exps, coeffs = convert_to_internal(representation, polynomials, ordering)

    # @info "wuw"
    # println(ring)
    # println(exps)
    # println(coeffs)
    #= check and set algorithm parameters =#
    metainfo = set_metaparameters(ring, ordering, certify, forsolve, linalg, rng)
    # now ring stores computation ordering
    # metainfo is now a struct to store target ordering

    iszerobasis = remove_zeros_from_input!(ring, exps, coeffs)
    iszerobasis && (return convert_to_output(ring, polynomials, exps, coeffs, metainfo))

    #= change input ordering if needed =#
    assure_ordering!(ring, exps, coeffs, metainfo)

    #= compute the groebner basis =#
    bexps, bcoeffs = groebner(ring, exps, coeffs, reduced, metainfo)

    # @info "uwow"
    # println(ring)
    # println(bexps)
    # println(bcoeffs)

    # ring contains ordering of computation, it is the requested ordering
    #= convert result back to representation of input =#
    convert_to_output(ring, polynomials, bexps, bcoeffs, metainfo)
end

#------------------------------------------------------------------------------
# Finite field groebner

#=
    Accepts arguments
    . ring - a polynomial ring where computation happens 
            (here, the characteristic should be nonzero)
=#
function groebner(
            ring::PolyRing,
            exps::Vector{Vector{M}},
            coeffs::Vector{Vector{C}},
            reduced::Bool,
            meta::GroebnerMetainfo{Rng}) where {M<:Monom, C<:CoeffFF, Rng<:Random.AbstractRNG}

    # select hashtable size
    tablesize = select_tablesize(ring, exps)

    basis, ht = initialize_structures(
                        ring, exps, coeffs, meta.rng, tablesize)

    f4!(ring, basis, ht, reduced, meta.linalg, meta.rng)

    gbexps = hash_to_exponents(basis, ht)
    gbexps, basis.coeffs
end

#------------------------------------------------------------------------------
# Rational field groebner

function groebner(
            ring::PolyRing,
            exps::Vector{Vector{M}},
            coeffs::Vector{Vector{C}},
            reduced::Bool,
            meta::GroebnerMetainfo) where {M<:Monom, C<:CoeffQQ}

    # we can mutate coeffs and exps here

    # select hashtable size
    tablesize = select_tablesize(ring, exps)
    @info "Selected tablesize $tablesize"

    # initialize hashtable and finite field generators structs
    gens_temp_ff, ht = initialize_structures_ff(ring, exps,
                                        coeffs, meta.rng, tablesize)
    gens_ff = deepcopy_basis(gens_temp_ff)

    # now hashtable is filled correctly,
    # and gens_temp_ff exponents are correct and in correct order.
    # gens_temp_ff coefficients are filled with random stuff and
    # gens_temp_ff.ch is 0

    # to store integer and rational coefficients of groebner basis
    coeffaccum  = CoeffAccum{BigInt, Rational{BigInt}}()
    # to store BigInt buffers and reduce overall memory usage
    coeffbuffer = CoeffBuffer()

    # scale coefficients of input to integers
    coeffs_zz = scale_denominators(coeffbuffer, coeffs)

    # keeps track of used prime numbers
    primetracker = PrimeTracker(coeffs_zz)

    #=
        exps, coeffs and coeffs_zz must be *not changed* during whole computation
    =#

    i = 1

    # copy basis so that we initial exponents dont get lost
    gens_ff = deepcopy_basis(gens_temp_ff)

    prime = nextluckyprime!(primetracker)
    @info "$i: selected lucky prime $prime"

    # perform reduction and store result in gens_ff
    reduce_modulo!(coeffbuffer, coeffs_zz, gens_ff.coeffs, prime)

    # do some things to ensure generators are correct
    cleanup_basis!(ring, gens_ff, prime)

    # compute groebner basis in finite field
    #=
    Need to make sure input invariants in f4! are satisfied, f4.jl for details
    =#

    tracer = Tracer()

    pairset = initialize_pairset(powertype(M))

    f4!(ring, gens_ff, tracer, pairset, ht, reduced, meta.linalg, meta.rng)
    # f4!(ring, gens_ff, ht, reduced, meta.linalg)

    # reconstruct into integers
    @info "CRT modulo ($(primetracker.modulo), $(prime))"

    # ADDED
    # global RECTIME
    # ADDED
    reconstruct_crt!(coeffbuffer, coeffaccum, primetracker, gens_ff.coeffs, prime)
    # reconstruct_crt!(coeffbuffer, coeffaccum, primetracker, gens_ff.coeffs, prime)

    # reconstruct into rationals
    @info "Reconstructing modulo $(primetracker.modulo)"
    # ADDED
    reconstruct_modulo!(coeffbuffer, coeffaccum, primetracker)
    # reconstruct_modulo!(coeffbuffer, coeffaccum, primetracker)

    correct = false
    # ADDED
    # global CORRTIME
    t = time()
    if correctness_check!(coeffaccum, coeffbuffer, primetracker, meta,
                            ring, exps, coeffs, coeffs_zz, gens_temp_ff, gens_ff, ht)
        @info "Success!"
        # CORRTIME += time() - t
        correct = true
    end
    # CORRTIME += time() - t

    gap = 1
    primegaps = (1, 1, 1, 1, 1)
    primemult = 2
    # if first prime was not successfull
    while !correct

        if i <= length(primegaps)
            gap = primegaps[i]
        else
            gap = gap*primemult
        end

        for j in 1:gap
            i += 1

            # copy basis so that initial exponents dont get lost
            gens_ff = deepcopy_basis(gens_temp_ff)
            prime = nextluckyprime!(primetracker)
            @info "$i: selected lucky prime $prime"
            # perform reduction and store result in gens_ff
            reduce_modulo!(coeffbuffer, coeffs_zz, gens_ff.coeffs, prime)
            # do some things to ensure generators are correct
            cleanup_basis!(ring, gens_ff, prime)
            # compute groebner basis in finite field
            #=
            Need to make sure input invariants in f4! are satisfied, f4.jl for details
            =#
            # ADDED
            # global F4TIME
            # ADDED

            # @error "tracer" tracer

            f4!(ring, gens_ff, tracer, pairset, ht, reduced, meta.linalg, meta.rng)
            # f4!(ring, gens_ff, ht, reduced, meta.linalg)
            # reconstruct into integers
            @info "CRT modulo ($(primetracker.modulo), $(prime))"
            # ADDED
            # global RECTIME
            # ADDED

            if meta.linalg === :prob
                # check if probabilistic linear algebra failed
                shape_control(coeffaccum.gb_coeffs_qq, gens_ff.coeffs)
            end

            # ADDED
            # prevcoeffs = deepcopy(coeffaccum.gb_coeffs_zz)

            reconstruct_crt!(coeffbuffer, coeffaccum, primetracker, gens_ff.coeffs, prime)
            # reconstruct_crt!(coeffbuffer, coeffaccum, primetracker, gens_ff.coeffs, prime)

            #=
            -----------------------------------
            Given r1, r2, and m1, m2, so that
                r1 ≡ 1/b1 (mod m1)
                r2 ≡ 1/b2 (mod m2)
            Determine if b1 = b2
            -----------------------------------

            Easy:
            b1 is inverse of r1 mod m1;
            b2 is inverse of r2 mod m2;

            then check equality
            =#

            #=
            -----------------------------------
            We go further.
            Given r1, r2, and m1, m2, so that
                r1 ≡ a1/b1 (mod m1)
                r2 ≡ a2/b2 (mod m2)
            Determine if a1 = a2 and b1 = b2
            -----------------------------------
            r1*b1 ≡ a1 (mod m1),
            r2*b2 ≡ a2 (mod m2),

            r1*b1 * r2*b2 = a1*a2 (mod m1*m2)
            =#

            #=
            if j != gap
                if crt_issame_check!(coeffbuffer, coeffaccum,
                                gens_ff.coeffs, primetracker.modulo, prime)
                    @info "CRT same check passed!"
                    break
                end
            end
            =#
        end

        # reconstruct into rationals
        @info "Reconstructing modulo $(primetracker.modulo)"
        # ADDED
        success = reconstruct_modulo!(coeffbuffer, coeffaccum, primetracker)
        if !success
            continue
        end

        # reconstruct_modulo!(coeffbuffer, coeffaccum, primetracker)
        # ADDED
        # global CORRTIME
        t = time()
        if correctness_check!(coeffaccum, coeffbuffer, primetracker, meta,
                                ring, exps, coeffs, coeffs_zz, gens_temp_ff, gens_ff, ht)
            @info "Success!"
            # CORRTIME += time() - t
            correct = true
        end
        # CORRTIME += time() - t
        #=
        if correctness_check!(coeffaccum, coeffbuffer, primetracker, meta,
                                ring, exps, coeffs, coeffs_zz, gens_temp_ff, gens_ff, ht)
            @info "Success!"
            break
        end
        =#

        if meta.linalg === :prob
            if i > 1024
                coeff_control(primetracker, coeffaccum)
            end
        end
    end

    #=
    if meta.usefglm
        fglm_f4!(ring, gens_ff, ht)
    end
    =#

    # ADDED
    # global PRIMES
    # PRIMES = i

    # normalize_coeffs!(gbcoeffs_qq)
    gb_exps = hash_to_exponents(gens_ff, ht)
    gb_exps, coeffaccum.gb_coeffs_qq
end

function assure_ordering!(ring, exps, coeffs, metainfo)
    if ring.ord != metainfo.computeord
        sort_input_to_change_ordering(exps, coeffs, metainfo.computeord)
    end
    ring.ord = metainfo.computeord
end

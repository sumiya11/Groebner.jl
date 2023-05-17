 
function groebner(
        polynomials::AbstractVector, 
        representation::RepresentationStyle, 
        reduced::Bool, 
        ordering::AbstractMonomialOrdering, 
        certify::Bool, 
        linalg::Symbol, 
        rng,
        maxpairs::Int)
    #= extract ring information, exponents and coefficients
       from input polynomials =#
    # Copies input, so that polynomials would not be changed itself.
    ring, exps, coeffs = convert_to_internal(representation, polynomials, ordering)

    #= check and set algorithm parameters =#
    metainfo = set_metaparameters(ring, ordering, certify, linalg, rng, maxpairs)
    # now ring stores computation ordering
    # metainfo is now a struct to store target ordering

    iszerobasis = remove_zeros_from_input!(ring, exps, coeffs)
    iszerobasis && (return convert_to_output(ring, polynomials, exps, coeffs, metainfo))

    #= change input ordering if needed =#
    newring = assure_ordering!(ring, exps, coeffs, metainfo.targetord)

    #= compute the groebner basis =#
    bexps, bcoeffs = groebner(newring, exps, coeffs, reduced, metainfo, maxpairs)

    # ring contains ordering of computation, it is the requested ordering
    #= convert result back to representation of input =#
    convert_to_output(newring, polynomials, bexps, bcoeffs, metainfo)
end

#------------------------------------------------------------------------------
# Finite field groebner

# groebner over integers modulo a prime is simple:
# just initialize and call f4 modulo a prime.

function groebner(
        ring::PolyRing,
        exps::Vector{Vector{M}},
        coeffs::Vector{Vector{C}},
        reduced::Bool,
        meta::GroebnerMetainfo{Rng},
        maxpairs) where {M<:Monom, C<:CoeffFF, Rng<:Random.AbstractRNG}

    # select hashtable size
    tablesize = select_tablesize(ring, exps)

    # initialize basis and hashtable structures
    basis, ht = initialize_structures(
                        ring, exps, coeffs, meta.rng, tablesize)

    f4!(ring, basis, ht, reduced, meta.linalg, meta.rng, maxpairs=maxpairs)

    # extract exponents from hashtable
    gbexps = hash_to_exponents(basis, ht)

    gbexps, basis.coeffs
end

#------------------------------------------------------------------------------
# Rational numbers groebner

# groebner over rationals is implemented roughly in the following way:
#
# k = 1
# while !(correctly reconstructed)
#   k = k*2
#   select a batch of small prime numbers p1..pk
#   compute a batch of finite field groebner bases gb1..gbk
#   reconstruct gb1..gbk to gb_zz modulo prod(p1..pk) with CRT
#   reconstruct gb_zz to gb_qq with rational reconstruction
#   if the basis gb_qq is correct, then break
# end
# return gb_qq
#

initial_batchsize() = 1
initial_gaps() = (1, 1, 1, 1, 1)
batchsize_multiplier() = 2

function groebner(
            ring::PolyRing,
            exps::Vector{Vector{M}},
            coeffs::Vector{Vector{C}},
            reduced::Bool,
            meta::GroebnerMetainfo,
            maxpairs) where {M<:Monom, C<:CoeffQQ}

    # we can mutate coeffs and exps here

    # select hashtable size
    tablesize = select_tablesize(ring, exps)
    @info "Selected tablesize $tablesize"

    # initialize hashtable and finite field basis structs
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

    # exps, coeffs and coeffs_zz **must be not changed** during whole computation
    i = 1

    # copy basis so that we initial exponents dont get lost
    gens_ff = deepcopy_basis(gens_temp_ff)

    prime = nextluckyprime!(primetracker)
    @info "$i: selected lucky prime $prime"

    # perform reduction and store result in gens_ff
    reduce_modulo!(coeffbuffer, coeffs_zz, gens_ff.coeffs, prime)

    # do some things to ensure generators are correct
    cleanup_basis!(ring, gens_ff, prime)

    # compute groebner basis in finite field.
    # Need to make sure input invariants in f4! are satisfied, see f4/f4.jl for details

    tracer = Tracer()
    pairset = initialize_pairset(powertype(M))

    f4!(ring, gens_ff, tracer, pairset, ht, reduced, meta.linalg, meta.rng, maxpairs)

    # reconstruct into integers
    @info "CRT modulo ($(primetracker.modulo), $(prime))"

    reconstruct_crt!(coeffbuffer, coeffaccum, primetracker, gens_ff.coeffs, prime)

    # reconstruct into rationals
    @info "Reconstructing modulo $(primetracker.modulo)"
    reconstruct_modulo!(coeffbuffer, coeffaccum, primetracker)

    correct = false
    if correctness_check!(coeffaccum, coeffbuffer, primetracker, meta,
                            ring, exps, coeffs, coeffs_zz, gens_temp_ff, gens_ff, ht)
        @info "Reconstructed successfuly"
        correct = true
    end
    
    # initial batch size
    batchsize = initial_batchsize()
    gaps = initial_gaps()
    # multiplier for the size of the batch of prime numbers
    multiplier = batchsize_multiplier()
    # if first prime was not successfull
    while !correct
        if i <= length(gaps)
            batchsize = gaps[i]
        else
            batchsize = batchsize*multiplier
        end

        for j in 1:batchsize
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
            # Need to make sure input invariants in f4! are satisfied, see f4/f4.jl for details

            f4!(ring, gens_ff, tracer, pairset, ht, reduced, meta.linalg, meta.rng, maxpairs)
            # reconstruct to integers
            @info "CRT modulo ($(primetracker.modulo), $(prime))"

            batchsize > 1 && basis_shape_control(coeffaccum.gb_coeffs_qq, gens_ff.coeffs)

            reconstruct_crt!(coeffbuffer, coeffaccum, primetracker, gens_ff.coeffs, prime)
        end

        # reconstruct to rationals
        @info "Reconstructing modulo $(primetracker.modulo)"
        success = reconstruct_modulo!(coeffbuffer, coeffaccum, primetracker)
        !success && continue

        if correctness_check!(coeffaccum, coeffbuffer, primetracker, meta,
                                ring, exps, coeffs, coeffs_zz, gens_temp_ff, gens_ff, ht)
            @info "Success!"
            correct = true
        end

        # i > 2^10 && basis_coeff_control(primetracker, coeffaccum)
    end

    #=
    # TODO
    if meta.usefglm
        fglm_f4!(ring, gens_ff, ht)
    end
    =#

    gb_exps = hash_to_exponents(gens_ff, ht)
    gb_exps, coeffaccum.gb_coeffs_qq
end


# Backend for `groebner`

# Proxy function dedicated to handling errors.
function _groebner(polynomials, kws::KeywordsHandler)
    # We try to guess efficient internal polynomial representation, i.e., a
    # suitable representation of monomials and coefficients.
    polynomial_repr = select_polynomial_representation(polynomials, kws)
    try
        # The backend is wrapped in a try/catch to catch exceptions that one can
        # hope to recover from (and, perhaps, restart the computation with safer
        # parameters).
        return _groebner(polynomials, kws, polynomial_repr)
    catch err
        if isa(err, ExponentVectorOverflow)
            # Exponent vector overflow might have happened. Restart the
            # computation with at least 32 bits per exponent.
            @log level = 0 "Possible overflow of exponent vector detected. Restarting with $(32) bits per exponent."
            polynomial_repr =
                select_polynomial_representation(polynomials, kws, hint=:large_exponents)
            return _groebner(polynomials, kws, polynomial_repr)
        else
            # Something bad happened.
            rethrow(err)
        end
    end
end

function _groebner(
    polynomials,
    kws::KeywordsHandler,
    representation::PolynomialRepresentation
)
    # Extract ring information, exponents, and coefficients from the input polynomials.
    # Convert these to an internal polynomial representation.
    # This must copy the input, so that `polynomials` itself is not modified.
    ring, monoms, coeffs = convert_to_internal(representation, polynomials, kws)
    # Fast path for the input of zeros
    # allzeros = remove_zeros_from_input!(ring, monoms, coeffs)
    # allzeros && return convert_to_output(ring, polynomials, monoms, coeffs, params)
    # Check and set parameters
    params = AlgorithmParameters(ring, kws)
    @log level = 1 "Selected parameters:\n$(params)"
    # NOTE: at this point, we already know the computation method we are going to use,
    # and the parameters are set.
    change_ordering_if_needed!(ring, monoms, coeffs, params)
    # Compute a groebner basis!
    gbmonoms, gbcoeffs = _groebner(ring, monoms, coeffs, params)
    # Convert result back to the representation of input
    convert_to_output(ring, polynomials, gbmonoms, gbcoeffs, params)
end

# Groebner basis over Z_p.
# Just calls f4 directly.
function _groebner(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffFF}
    # NOTE: we can mutate ring, monoms, and coeffs here.
    @log level = 1 "Backend: F4 over Z_$(ring.ch)"
    basis, pairset, hashtable = initialize_structs(ring, monoms, coeffs, params)
    tracer = Tracer()
    f4!(ring, basis, pairset, hashtable, tracer, params)
    # Extract monomials and coefficients from basis and hashtable
    gbmonoms, gbcoeffs = extract_monoms_coeffs(basis, hashtable)
    gbmonoms, gbcoeffs
end

# Groebner basis over Q.
# GB over the rationals uses modular computation.
function _groebner(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    # NOTE: we can mutate ring, monoms, and coeffs here.
    if params.strategy === :learn_and_apply
        _groebner_learn_and_apply(ring, monoms, coeffs, params)
    else
        @assert params.strategy === :classic_modular
        _groebner_classic_modular(ring, monoms, coeffs, params)
    end
end

function _groebner_learn_and_apply(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    @log level = 1 "Backend: multi-modular tracing (learn-apply)"
    monoms, coeffs
end

function _groebner_classic_modular(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    @log level = 1 "Backend: classic multi-modular"
    basis, pairset, hashtable, tracer = initialize_structs(ring, monoms, coeffs, params)
    gens_ff = deepcopy_basis(gens_temp_ff)

    # now hashtable is filled correctly,
    # and gens_temp_ff exponents are correct and in correct order.
    # gens_temp_ff coefficients are filled with random stuff and
    # gens_temp_ff.ch is 0

    # to store integer and rational coefficients of groebner basis
    coeffaccum = CoeffAccum{BigInt, Rational{BigInt}}()
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
    # @info "$i: selected lucky prime $prime"

    # perform reduction and store result in gens_ff
    reduce_modulo!(coeffbuffer, coeffs_zz, gens_ff.coeffs, prime)

    # do some things to ensure generators are correct
    cleanup_basis!(ring, gens_ff, prime)

    # compute groebner basis in finite field.
    # Need to make sure input invariants in f4! are satisfied, see f4/f4.jl for details

    tracer = Tracer()
    pairset = initialize_pairset(powertype(M))

    f4!(ring, gens_ff, tracer, pairset, ht, params)

    # reconstruct into integers
    # @info "CRT modulo ($(primetracker.modulo), $(prime))"

    reconstruct_crt!(coeffbuffer, coeffaccum, primetracker, gens_ff.coeffs, prime)

    # reconstruct into rationals
    # @info "Reconstructing modulo $(primetracker.modulo)"
    reconstruct_modulo!(coeffbuffer, coeffaccum, primetracker)

    correct = false
    if correctness_check!(
        coeffaccum,
        coeffbuffer,
        primetracker,
        meta,
        ring,
        exps,
        coeffs,
        coeffs_zz,
        gens_temp_ff,
        gens_ff,
        ht
    )
        # @info "Reconstructed successfuly"
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
            batchsize = batchsize * multiplier
        end

        for j in 1:batchsize
            i += 1

            # copy basis so that initial exponents dont get lost
            gens_ff = deepcopy_basis(gens_temp_ff)
            prime = nextluckyprime!(primetracker)
            # @info "$i: selected lucky prime $prime"
            # perform reduction and store result in gens_ff
            reduce_modulo!(coeffbuffer, coeffs_zz, gens_ff.coeffs, prime)
            # do some things to ensure generators are correct
            cleanup_basis!(ring, gens_ff, prime)
            # compute groebner basis in finite field
            # Need to make sure input invariants in f4! are satisfied, see f4/f4.jl for details

            f4!(
                ring,
                gens_ff,
                tracer,
                pairset,
                ht,
                reduced,
                meta.linalg,
                meta.rng,
                maxpairs
            )
            # reconstruct to integers
            # @info "CRT modulo ($(primetracker.modulo), $(prime))"

            batchsize > 1 && basis_shape_control(coeffaccum.gb_coeffs_qq, gens_ff.coeffs)

            reconstruct_crt!(coeffbuffer, coeffaccum, primetracker, gens_ff.coeffs, prime)
        end

        # reconstruct to rationals
        # @info "Reconstructing modulo $(primetracker.modulo)"
        success = reconstruct_modulo!(coeffbuffer, coeffaccum, primetracker)
        !success && continue

        if correctness_check!(
            coeffaccum,
            coeffbuffer,
            primetracker,
            meta,
            ring,
            exps,
            coeffs,
            coeffs_zz,
            gens_temp_ff,
            gens_ff,
            ht
        )
            # @info "Success!"
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

# This file is a part of Groebner.jl. License is GNU GPL v2.

### 
# Backend for `groebner_with_change_matrix`

# Proxy function for handling exceptions.
function _groebner_with_change_matrix0(polynomials, kws::KeywordArguments)
    # We try to select an efficient internal polynomial representation, i.e., a
    # suitable representation of monomials and coefficients.
    polynomial_repr = io_select_polynomial_representation(polynomials, kws)
    try
        # The backend is wrapped in a try/catch to catch exceptions that one can
        # hope to recover from (and, perhaps, restart the computation with safer
        # parameters).
        return _groebner_with_change_matrix1(polynomials, kws, polynomial_repr)
    catch err
        if isa(err, MonomialDegreeOverflow)
            @log :info """
            Possible overflow of exponent vector detected. 
            Restarting with at least $(32) bits per exponent."""
            polynomial_repr =
                io_select_polynomial_representation(polynomials, kws, hint=:large_exponents)
            return _groebner_with_change_matrix1(polynomials, kws, polynomial_repr)
        else
            # Something bad happened.
            rethrow(err)
        end
    end
end

function _groebner_with_change_matrix1(polynomials, kws::KeywordArguments, representation)
    # Extract ring information, exponents, and coefficients from input
    # polynomials. Convert these to an internal polynomial representation. 
    # NOTE: This must copy the input, so that input `polynomials` is never
    # modified.
    # NOTE: The body of this function is type-unstable (by design)
    ring, var_to_index, monoms, coeffs =
        io_convert_to_internal(representation, polynomials, kws)

    # Check and set parameters and monomial ordering
    params = AlgorithmParameters(ring, representation, kws)
    ring, _ = io_set_monomial_ordering!(ring, var_to_index, monoms, coeffs, params)

    # Fast path for the input of zeros
    if isempty(monoms)
        @log :misc "Input consisting of zero polynomials. Returning zero."
        changematrix = io_convert_changematrix_to_output(
            ring,
            polynomials,
            0,
            [empty(monoms)],
            [empty(coeffs)],
            params
        )
        return io_convert_to_output(ring, polynomials, monoms, coeffs, params), changematrix
    end

    if params.homogenize
        # this also performs saturation w.r.t. the homogenizing variable
        _, ring, monoms, coeffs = homogenize_generators!(ring, monoms, coeffs, params)
    end

    # Compute a groebner basis!
    gbmonoms, gbcoeffs, matrix_monoms, matrix_coeffs =
        _groebner_with_change_matrix2(ring, monoms, coeffs, params)

    if params.homogenize
        ring, gbmonoms, gbcoeffs =
            dehomogenize_generators!(ring, gbmonoms, gbcoeffs, params)
    end

    # Convert result back to the representation of input
    basis = io_convert_to_output(ring, polynomials, gbmonoms, gbcoeffs, params)
    changematrix = io_convert_changematrix_to_output(
        ring,
        polynomials,
        length(monoms),
        matrix_monoms,
        matrix_coeffs,
        params
    )

    basis, changematrix
end

###
# Groebner basis over Z_p.
# Just calls f4 directly.

@timeit function _groebner_with_change_matrix2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffZp}
    # NOTE: we can mutate ring, monoms, and coeffs here.
    @log :misc "Backend: F4 over Z_$(ring.ch) (with change matrix)"
    # NOTE: the sorting of input polynomials is not deterministic across
    # different Julia versions when sorting only w.r.t. the leading term
    basis, pairset, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    tracer = TinyTraceF4()
    f4!(ring, basis, pairset, hashtable, tracer, params)
    # Extract monomials and coefficients from basis and hashtable
    gbmonoms, gbcoeffs = basis_export_data(basis, hashtable)
    matrix_monoms, matrix_coeffs =
        basis_changematrix_export(basis, hashtable, length(monoms))
    gbmonoms, gbcoeffs, matrix_monoms, matrix_coeffs
end

###
# Groebner basis over Q.

# GB over the rationals uses modular computation.
@timeit function _groebner_with_change_matrix2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    _groebner_with_change_classic_modular(ring, monoms, coeffs, params)
end

###
# Classic multi-modular strategy

function _groebner_with_change_classic_modular(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    # NOTE: we can mutate ring, monoms, and coeffs here.
    @log :misc "Backend: classic multi-modular F4 (with change matrix)"

    # Initialize supporting structs
    state = GroebnerState{BigInt, C, CoeffModular}(params)
    # Initialize F4 structs
    basis, pairset, hashtable =
        f4_initialize_structs(ring, monoms, coeffs, params, normalize_input=false)
    tracer = TinyTraceF4(params)
    crt_algorithm = params.crt_algorithm

    # Scale the input coefficients to integers to speed up the subsequent search
    # for lucky primes
    @log :debug "Input polynomials" basis
    @log :misc "Clearing the denominators of the input polynomials"
    basis_zz = clear_denominators!(state.buffer, basis, deepcopy=false)
    @log :debug "Integer coefficients are" basis_zz.coeffs

    # Handler for lucky primes
    luckyprimes = LuckyPrimes(basis_zz.coeffs)
    prime = next_lucky_prime!(luckyprimes)
    @log :misc "The first lucky prime is $prime"
    @log :misc "Reducing input generators modulo $prime"

    # Perform reduction modulo prime and store result in basis_ff
    ring_ff, basis_ff = reduce_modulo_p!(state.buffer, ring, basis_zz, prime, deepcopy=true)
    @log :debug "Reduced coefficients are" basis_ff.coeffs

    @log :debug "Before F4" basis_ff
    params_zp = params_mod_p(params, prime)
    f4!(ring_ff, basis_ff, pairset, hashtable, tracer, params_zp)
    @log :debug "After F4:" basis_ff
    # NOTE: basis_ff may not own its coefficients, one should not mutate it
    # directly further in the code

    push!(state.gb_coeffs_ff_all, basis_ff.coeffs)

    changematrix_monoms, changematrix_coeffs =
        basis_changematrix_export(basis_ff, hashtable, length(monoms))
    push!(state.changematrix_coeffs_ff_all, changematrix_coeffs)

    # Reconstruct coefficients and write results to the accumulator.
    # CRT reconstrction is trivial here.
    @log :misc "Reconstructing coefficients from Z_$prime to QQ"
    full_simultaneous_crt_reconstruct!(state, luckyprimes)
    full_simultaneous_crt_reconstruct_changematrix!(state, luckyprimes)

    success_reconstruct = full_rational_reconstruct!(state, luckyprimes, params.use_flint)
    @log :debug "Reconstructed coefficients" state.gb_coeffs_qq
    @log :misc "Successfull reconstruction: $success_reconstruct"

    changematrix_success_reconstruct =
        full_rational_reconstruct_changematrix!(state, luckyprimes, params.use_flint)
    @log :misc "Successfull change matrix reconstruction: $success_reconstruct"

    correct_basis = false
    if success_reconstruct && changematrix_success_reconstruct
        @log :misc "Verifying the correctness of reconstruction"
        correct_basis = correctness_check!(
            state,
            luckyprimes,
            ring_ff,
            basis,
            basis_zz,
            basis_ff,
            hashtable,
            params
        )
        @log :misc "Passed correctness check: $correct_basis"
        # At this point, if the constructed basis is correct, we return it.
        if correct_basis
            # take monomials from the basis modulo a prime
            gb_monoms, _ = basis_export_data(basis_ff, hashtable)
            # take coefficients from the reconstrcted basis
            gb_coeffs_qq = state.gb_coeffs_qq
            changematrix_coeffs_qq = state.changematrix_coeffs_qq
            return gb_monoms, gb_coeffs_qq, changematrix_monoms, changematrix_coeffs_qq
        end
    end

    # At this point, either the reconstruction or the correctness check failed.
    # Continue to compute Groebner bases modulo different primes in batches. 
    primes_used = 1
    batchsize = 1
    batchsize_scaling = 0.10
    @log :misc """
    Preparing to compute bases in batches.. 
    The initial size of the batch is $batchsize. 
    The batch scale factor is $batchsize_scaling."""

    if !tracer.ready_to_use
        @log :misc """
          The tracer is disabled until the shape of the basis is not determined via majority vote.
          The threshold for the majority vote is $(params.majority_threshold)
          """
    end

    # CRT and rational reconstrction settings
    indices_selection = Vector{Tuple{Int, Int}}(undef, length(state.gb_coeffs_zz))
    k = 1
    for i in 1:length(state.gb_coeffs_zz)
        l = length(state.gb_coeffs_zz[i])
        nl = max(isqrt(l) - 1, 1)
        if isone(l)
            continue
        end
        for j in 1:nl
            if k > length(indices_selection)
                resize!(indices_selection, 2 * length(indices_selection))
            end
            indices_selection[k] = (i, rand(2:l))
            k += 1
        end
    end
    resize!(indices_selection, k - 1)
    unique!(indices_selection)

    partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)

    iters = 0
    while !correct_basis
        @log :misc "Iteration # $iters of modular Groebner"

        for j in 1:batchsize
            prime = next_lucky_prime!(luckyprimes)
            @log :debug "The lucky prime is $prime"
            @log :debug "Reducing input generators modulo $prime"

            # Perform reduction modulo prime and store result in basis_ff
            ring_ff, basis_ff =
                reduce_modulo_p!(state.buffer, ring, basis_zz, prime, deepcopy=true)
            params_zp = params_mod_p(params, prime)

            f4!(ring_ff, basis_ff, pairset, hashtable, tracer, params_zp)

            push!(state.gb_coeffs_ff_all, basis_ff.coeffs)
            _, changematrix_coeffs =
                basis_changematrix_export(basis_ff, hashtable, length(monoms))
            push!(state.changematrix_coeffs_ff_all, changematrix_coeffs)

            if !majority_vote!(state, basis_ff, tracer, params)
                @log :debug "Majority vote is not conclusive, aborting reconstruction!"
                continue
            end

            if crt_algorithm === :incremental
                partial_incremental_crt_reconstruct!(state, luckyprimes, indices_selection)
            end
            primes_used += 1
        end
        if crt_algorithm === :simultaneous
            partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)
        end

        @log :misc "Partially reconstructing coefficients to QQ"
        @log :debug "Partially reconstructing coefficients from Z_$(luckyprimes.modulo * prime) to QQ"
        success_reconstruct = partial_rational_reconstruct!(
            state,
            luckyprimes,
            indices_selection,
            params.use_flint
        )
        @log :misc "Partial reconstruction successfull: $success_reconstruct"
        @log :misc """
          Used $(length(luckyprimes.primes)) primes in total over $(iters + 1) iterations.
          The current batch size is $batchsize.
          """

        if !success_reconstruct
            @log :misc "Partial rational reconstruction failed"
            iters += 1
            batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
            continue
        end

        if params.heuristic_check
            success_check =
                heuristic_correctness_check(state.selected_coeffs_qq, luckyprimes.modulo)
            if !success_check
                @log :misc "Heuristic check failed for partial reconstruction"
                iters += 1
                batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
                continue
            end
        end

        # Perform full reconstruction
        @log :misc "Performing full CRT and rational reconstruction.."
        if crt_algorithm === :incremental
            full_incremental_crt_reconstruct!(state, luckyprimes)
        else
            full_simultaneous_crt_reconstruct!(state, luckyprimes)
        end
        success_reconstruct =
            full_rational_reconstruct!(state, luckyprimes, params.use_flint)

        if !success_reconstruct
            @log :misc "Full reconstruction failed"
            iters += 1
            batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
            continue
        end

        full_simultaneous_crt_reconstruct_changematrix!(state, luckyprimes)
        success_reconstruct_changematrix =
            full_rational_reconstruct_changematrix!(state, luckyprimes, params.use_flint)
        if !success_reconstruct_changematrix
            @log :misc "Failed to reconstruct the change matrix"
            iters += 1
            batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
            continue
        end

        correct_basis = correctness_check!(
            state,
            luckyprimes,
            ring_ff,
            basis,
            basis_zz,
            basis_ff,
            hashtable,
            params
        )

        iters += 1
        batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
    end

    @log :misc "Correctness check passed!"
    @log :misc "Used $(length(luckyprimes.primes)) primes in total over $(iters) iterations"

    # Construct the output basis.
    # Take monomials from the basis modulo a prime
    gb_monoms, _ = basis_export_data(basis_ff, hashtable)
    # Take coefficients from the reconstructed basis
    gb_coeffs_qq = state.gb_coeffs_qq
    changematrix_coeffs_qq = state.changematrix_coeffs_qq

    return gb_monoms, gb_coeffs_qq, changematrix_monoms, changematrix_coeffs_qq
end

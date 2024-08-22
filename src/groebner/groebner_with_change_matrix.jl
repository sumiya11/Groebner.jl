# This file is a part of Groebner.jl. License is GNU GPL v2.

### 
# Backend for `groebner_with_change_matrix`

function groebner_with_change_matrix0(polynomials, options)
    isempty(polynomials) && throw(DomainError("Empty input."))
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    gb_monoms, gb_coeffs, matrix_monoms, matrix_coeffs =
        _groebner_with_change_matrix1(ring, monoms, coeffs, options)
    result = io_convert_ir_to_polynomials(ring, polynomials, gb_monoms, gb_coeffs, options)
    rows = [
        io_convert_ir_to_polynomials(ring, polynomials, m, c, options) for
        (m, c) in zip(matrix_monoms, matrix_coeffs)
    ]
    matrix = collect(transpose(reduce(hcat, rows)))
    result, matrix
end

function groebner_with_change_matrix1(ring, monoms, coeffs, options)
    ring, monoms, coeffs = ir_ensure_assumptions(ring, monoms, coeffs)
    _groebner_with_change_matrix1(ring, monoms, coeffs, options)
end

function _groebner_with_change_matrix1(ring::PolyRing, monoms, coeffs, options)
    try
        params = AlgorithmParameters(ring, options)
        return __groebner_with_change_matrix1(ring, monoms, coeffs, params)
    catch err
        if isa(err, MonomialDegreeOverflow)
            @log :info """
            Possible overflow of exponent vector detected. 
            Restarting with at least 32 bits per exponent.""" maxlog = 1
            params = AlgorithmParameters(ring, options; hint=:large_exponents)
            return __groebner_with_change_matrix1(ring, monoms, coeffs, params)
        else
            # Something bad happened.
            rethrow(err)
        end
    end
end

function __groebner_with_change_matrix1(ring, monoms, coeffs, params)
    @invariant ir_is_valid(ring, monoms, coeffs)
    _, ring2, monoms2, coeffs2 =
        io_convert_ir_to_internal(ring, monoms, coeffs, params, params.representation)
    gb_monoms2, gb_coeffs2, matrix_monoms2, matrix_coeffs2 =
        groebner_with_change_matrix2(ring2, monoms2, coeffs2, params)
    gb_monoms, gb_coeffs = io_convert_internal_to_ir(ring2, gb_monoms2, gb_coeffs2, params)
    matrix_monoms_and_coeffs = [
        io_convert_internal_to_ir(ring2, _matrix_monoms2, _matrix_coeffs2, params) for
        (_matrix_monoms2, _matrix_coeffs2) in zip(matrix_monoms2, matrix_coeffs2)
    ]
    matrix_monoms, matrix_coeffs =
        map(first, matrix_monoms_and_coeffs), map(last, matrix_monoms_and_coeffs)
    gb_monoms, gb_coeffs, matrix_monoms, matrix_coeffs
end

function groebner_with_change_matrix2(ring, monoms, coeffs, params)
    nonzero_indices = findall(!io_iszero_coeffs, coeffs)
    zero_indices = setdiff(collect(1:length(monoms)), nonzero_indices)

    __monoms = monoms
    __coeffs = coeffs
    _monoms = filter(!isempty, monoms)
    _coeffs = filter(!isempty, coeffs)
    if isempty(_monoms)
        return [monoms[1]],
        [coeffs[1]],
        [[monoms[1] for _ in 1:length(monoms)]],
        [[coeffs[1] for _ in 1:length(coeffs)]]
    end
    monoms, coeffs = _monoms, _coeffs

    if params.homogenize
        _, ring, monoms, coeffs = homogenize_generators!(ring, monoms, coeffs, params)
    end

    gb_monoms, gb_coeffs, matrix_monoms, matrix_coeffs =
        _groebner_with_change_matrix2(ring, monoms, coeffs, params)

    if params.homogenize
        ring, gb_monoms, gb_coeffs =
            dehomogenize_generators!(ring, gb_monoms, gb_coeffs, params)
    end

    if !isempty(zero_indices)
        for i in 1:length(matrix_monoms)
            _matrix_monoms_i = resize!(similar(matrix_monoms[i]), length(__monoms))
            _matrix_coeffs_i = resize!(similar(matrix_coeffs[i]), length(__monoms))
            zero_indices = setdiff(collect(1:length(__monoms)), nonzero_indices)
            _matrix_monoms_i .= [__monoms[zero_indices[1]]]
            _matrix_coeffs_i .= [__coeffs[zero_indices[1]]]
            for j in 1:length(nonzero_indices)
                _matrix_monoms_i[nonzero_indices[j]] = matrix_monoms[i][j]
                _matrix_coeffs_i[nonzero_indices[j]] = matrix_coeffs[i][j]
            end
            matrix_monoms[i] = _matrix_monoms_i
            matrix_coeffs[i] = _matrix_coeffs_i
        end
    end

    gb_monoms, gb_coeffs, matrix_monoms, matrix_coeffs
end

###
# Groebner basis over Z_p. Calls F4 directly.

function _groebner_with_change_matrix2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffZp}
    basis, pairset, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    f4!(ring, basis, pairset, hashtable, params)
    gbmonoms, gbcoeffs = basis_export_data(basis, hashtable)
    matrix_monoms, matrix_coeffs =
        basis_changematrix_export(basis, hashtable, length(monoms))
    gbmonoms, gbcoeffs, matrix_monoms, matrix_coeffs
end

###
# Groebner basis over Q.

# GB over the rationals uses modular computation.
function _groebner_with_change_matrix2(
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

    # Initialize supporting structs
    state = ModularState{BigInt, C, CoeffModular}(params)
    basis, pairset, hashtable =
        f4_initialize_structs(ring, monoms, coeffs, params, make_monic=false)

    basis_zz = clear_denominators!(basis, deepcopy=false)

    # Handler for lucky primes
    luckyprimes = LuckyPrimes(basis_zz.coeffs)
    prime = primes_next_lucky_prime!(luckyprimes)

    # Perform reduction modulo prime and store result in basis_ff
    ring_ff, basis_ff = reduce_modulo_p!(ring, basis_zz, prime, deepcopy=true)

    params_zp = params_mod_p(params, prime)
    f4!(ring_ff, basis_ff, pairset, hashtable, params_zp)
    # NOTE: basis_ff may not own its coefficients, one should not mutate it
    # directly further in the code

    push!(state.gb_coeffs_ff_all, basis_ff.coeffs)

    changematrix_monoms, changematrix_coeffs =
        basis_changematrix_export(basis_ff, hashtable, length(monoms))
    push!(state.changematrix_coeffs_ff_all, changematrix_coeffs)

    # Reconstruct coefficients and write results to the accumulator.
    # CRT reconstrction is trivial here.
    modular_crt_full!(state, luckyprimes)
    full_simultaneous_crt_reconstruct_changematrix!(state, luckyprimes)

    success_reconstruct = modular_ratrec_vec_full!(state, luckyprimes)

    changematrix_success_reconstruct =
        full_rational_reconstruct_changematrix!(state, luckyprimes)

    correct_basis = false
    if success_reconstruct && changematrix_success_reconstruct
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
        # At this point, if the constructed basis is correct, we return it.
        if correct_basis
            gb_monoms, _ = basis_export_data(basis_ff, hashtable)
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

    tot, witness_set = modular_witness_set(state.gb_coeffs_zz, params)

    modular_crt_partial!(state, luckyprimes, witness_set)

    iters = 0
    while !correct_basis
        for j in 1:batchsize
            prime = primes_next_lucky_prime!(luckyprimes)

            ring_ff, basis_ff = reduce_modulo_p!(ring, basis_zz, prime, deepcopy=true)
            params_zp = params_mod_p(params, prime)

            f4!(ring_ff, basis_ff, pairset, hashtable, params_zp)

            push!(state.gb_coeffs_ff_all, basis_ff.coeffs)
            _, changematrix_coeffs =
                basis_changematrix_export(basis_ff, hashtable, length(monoms))
            push!(state.changematrix_coeffs_ff_all, changematrix_coeffs)

            if !majority_vote!(state, basis_ff, params)
                continue
            end

            primes_used += 1
        end
        modular_crt_partial!(state, luckyprimes, witness_set)

        success_reconstruct =
            partial_rational_reconstruct!(state, luckyprimes, witness_set)

        if !success_reconstruct
            iters += 1
            batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
            continue
        end

        if params.heuristic_check
            success_check =
                heuristic_correctness_check(state.selected_coeffs_qq, luckyprimes.modulo)
            if !success_check
                iters += 1
                batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
                continue
            end
        end

        # Perform full reconstruction
        modular_crt_full!(state, luckyprimes)

        success_reconstruct = modular_ratrec_vec_full!(state, luckyprimes)

        if !success_reconstruct
            iters += 1
            batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
            continue
        end

        full_simultaneous_crt_reconstruct_changematrix!(state, luckyprimes)
        success_reconstruct_changematrix =
            full_rational_reconstruct_changematrix!(state, luckyprimes)
        if !success_reconstruct_changematrix
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

    # Construct the output basis.
    gb_monoms, _ = basis_export_data(basis_ff, hashtable)
    gb_coeffs_qq = state.gb_coeffs_qq
    changematrix_coeffs_qq = state.changematrix_coeffs_qq

    gb_monoms, gb_coeffs_qq, changematrix_monoms, changematrix_coeffs_qq
end

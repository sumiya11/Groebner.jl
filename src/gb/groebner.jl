
function groebner(
        polynomials::AbstractVector,
        representation::RepresentationStyle, 
        reduced::Bool, 
        ordering::Symbol, 
        threading::Bool,
        certify::Bool, 
        forsolve::Bool, 
        linalg::Symbol, 
        rng)
    #= extract ring information, exponents and coefficients
       from input polynomials =#
    # Copies input, so that polynomials would not be changed itself.
    ring, exps, coeffs = convert_to_internal(representation, polynomials, ordering)

    #= check and set algorithm parameters =#
    metainfo = set_metaparameters(ring, ordering, certify, forsolve, linalg, rng)
    # now ring stores computation ordering
    # metainfo is now a struct to store target ordering

    iszerobasis = remove_zeros_from_input!(ring, exps, coeffs)
    iszerobasis && (return convert_to_output(ring, polynomials, exps, coeffs, metainfo))

    #= change input ordering if needed =#
    assure_ordering!(ring, exps, coeffs, metainfo.targetord)

    #= compute the groebner basis =#
    bexps, bcoeffs = groebner(ring, exps, coeffs, reduced, threading, metainfo)
        
    # ring contains ordering of computation, it is the requested ordering
    #= convert result back to representation of input =#
    convert_to_output(ring, polynomials, bexps, bcoeffs, metainfo)
end

#------------------------------------------------------------------------------
# Finite field groebner

# groebner over integers modulo a prime is simple:
# just initialize and call f4 modulo a prime.

# Parallelism is not used for Groebner over a finite field

# Given the hashed monomials of the basis (in gens_ff),
# the hashtable that hashes the exponents (in ht),
# and the coefficients of the basis,
# Returns the pair (exponents, coefficients) of the basis
function construct_final_basis(gens_ff, ht, coeffs)
    gb_exps = hash_to_exponents(gens_ff, ht)
    gb_exps, coeffs
end

function groebner(
        ring::PolyRing,
        exps::Vector{Vector{M}},
        coeffs::Vector{Vector{C}},
        reduced::Bool,
        threading::Bool,
        meta::GroebnerMetainfo{Rng}) where {M<:Monom, C<:CoeffFF, Rng<:Random.AbstractRNG}

    # select hashtable size
    tablesize = select_tablesize(ring, exps)

    # initialize basis and hashtable structures
    basis, ht = initialize_structures(
                        ring, exps, coeffs, meta.rng, tablesize)

    f4!(ring, basis, ht, reduced, meta.linalg, meta.rng)

    # extract exponents from hashtable,
    # and return it together with coefficients
    construct_final_basis(basis, ht, basis.coeffs)
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

# The story with parallelism with Groebner over rationals.
# 
# 

#
#
function distribute_things_per_threads!(nbases_per_thread::Vector{Int}, batchsize::Int)
    n_active_threads = nthreads()
    q_per_thread = div(batchsize, n_active_threads)
    residual = 0
    if q_per_thread < 1
        n_active_threads = batchsize
        residual = batchsize
    elseif q_per_thread >= 1
        residual = batchsize - n_active_threads * q_per_thread
    end
    for i in 1:n_active_threads
        nbases_per_thread[i] = 0
        nbases_per_thread[i] = q_per_thread + (residual > 0)
        residual -= 1
    end
    n_active_threads
end

#
# 
function initial_batchsize(threading::Bool)
    if threading
        nthreads()
    else
        1
    end
end
#
#
function initial_gaps(threading::Bool)
    if threading
        ()
    else
        (1, 1, 1, 1, 1)
    end
end
#
#
function batchsize_multiplier(threading)
    if threading
        1.0
    else
        1.5
    end
end

#
#
function compute_single_basis(
        ring, gens_temp_ff, tracer, pairset,
        ht, reduced, meta, rng, 
        coeffs_zz, used_prime, 
        coeffbuffer)
    # copy basis so that initial exponents dont get lost
    gens_ff = deepcopy_basis(gens_temp_ff)
    # perform reduction and store result in gens_ff
    reduce_modulo!(coeffbuffer, coeffs_zz, gens_ff.coeffs, used_prime)
    # do some things to ensure generators are correct
    cleanup_basis!(ring, gens_ff, used_prime)
    # compute groebner basis in finite field
    # Need to make sure input invariants in f4! are satisfied, see f4/f4.jl for details
    f4!(ring, gens_ff, tracer, pairset, ht, reduced, meta.linalg, rng)
    gens_ff
end

#
#
function compute_bases_in_thread(
                ring, gens_temp_ff, tracer, pairset,
                ht, reduced, meta, rng, 
                coeffs_zz, used_primes, computed_bases, 
                coeffbuffer, from_idx, to_idx)
    for i in from_idx:to_idx
        gens_ff = deepcopy_basis(gens_temp_ff)
        prime = used_primes[i]
        reduce_modulo!(coeffbuffer, coeffs_zz, gens_ff.coeffs, prime)
        cleanup_basis!(ring, gens_ff, prime)
        f4!(ring, gens_ff, tracer, pairset, ht, reduced, meta.linalg, rng)
        computed_bases[i] = gens_ff
    end
    nothing
end

function groebner(
            ring::PolyRing,
            exps::Vector{Vector{M}},
            coeffs::Vector{Vector{C}},
            reduced::Bool,
            threading::Bool,
            meta::GroebnerMetainfo) where {M<:Monom, C<:CoeffQQ}
    # We can mutate coeffs and exps here.
    # Select hashtable size
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

    f4!(ring, gens_ff, tracer, pairset, ht, reduced, meta.linalg, meta.rng)

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
    # the total number of primes used
    nprimes_used = 1

    # if the first prime is successfull
    if correct
        return construct_final_basis(gens_ff, ht, coeffaccum.gb_coeffs_qq)
    end
    
    # initial batch size
    batchsize = initial_batchsize(threading)
    gaps = initial_gaps(threading)
    # multiplier for the size of the batch of prime numbers
    multiplier = batchsize_multiplier(threading)
    
    # the size of the next batch of prime numbers
    next_batchsize(i, sz) = i <= length(gaps) ? gaps[i] : floor(Int, sz*multiplier)

    if threading
        # @warn "Using $(nthreads()) threads."
        if nthreads() == 1
            @info "Groebner called with `threading=true`, 
                    but there is only one active thread."
        end

        # use Threads.resize_nthreads! ?
        nbases_per_thread = zeros(Int, nthreads())
        rings_per_thread = [deepcopy(ring) for _ in 1:nthreads()]
        coeff_buffers_per_thread = [CoeffBuffer() for _ in 1:nthreads()]
        pairset_per_thread = [initialize_pairset(powertype(M)) for _ in 1:nthreads()]
        hashtable_per_thread = [threadsafe_copy(ht) for _ in 1:nthreads()]
        rng_per_thread = [deepcopy(meta.rng) for _ in 1:nthreads()]
        
        futures = Vector{Task}(undef, nthreads())

        while !correct
            nprimes_used += batchsize
            batchsize = next_batchsize(nprimes_used, batchsize)

            # Distribute the bases over available threads evenly
            n_active_threads = distribute_things_per_threads!(nbases_per_thread, batchsize)

            computed_bases = Vector{typeof(gens_temp_ff)}(undef, batchsize)
            used_primes = Vector{CoeffModular}(undef, batchsize)
            for j in 1:batchsize
                prime = nextluckyprime!(primetracker)
                used_primes[j] = prime
            end

            # Spawn n_active_threads threads,
            # each thread will compute bases from from_idx to to_idx
            # (to_idx - from_idx + 1 bases in total)
            absolute_idx = 1
            for th in 1:n_active_threads
                from_idx = absolute_idx
                to_idx = from_idx + nbases_per_thread[th] - 1
                absolute_idx = to_idx + 1
                futures[th] = @tspawnat th compute_bases_in_thread(
                    rings_per_thread[th],
                    gens_temp_ff,
                    tracer,
                    pairset_per_thread[th],
                    hashtable_per_thread[th],
                    reduced,
                    meta, 
                    rng_per_thread[th],
                    coeffs_zz,
                    used_primes,
                    computed_bases,
                    coeff_buffers_per_thread[th],
                    from_idx, to_idx
                )  
            end

            for j in 1:n_active_threads
                wait(futures[j])
            end

            npolys = length(gens_ff.coeffs)
            # Distribute the polynomials over available threads evenly
            n_active_threads = distribute_things_per_threads!(nbases_per_thread, npolys)
            threaded_reconstruct_coeffs!(
                nbases_per_thread, n_active_threads,
                coeff_buffers_per_thread, coeffaccum, primetracker, 
                used_primes, computed_bases
            )

            # for j in 1:batchsize
            #     prime = used_primes[j]
            #     gens_ff = computed_bases[j]
            #     batchsize > 1 && basis_shape_control(coeffaccum.gb_coeffs_qq, gens_ff.coeffs)
            #     reconstruct_crt!(coeffbuffer, coeffaccum, primetracker, gens_ff.coeffs, prime)
            # end
            
            # @info "Reconstructing modulo $(primetracker.modulo)"
            # success = reconstruct_modulo!(coeffbuffer, coeffaccum, primetracker)
            # !success && continue

            if correctness_check!(coeffaccum, coeffbuffer, primetracker, meta,
                                    ring, exps, coeffs, coeffs_zz, gens_temp_ff, gens_ff, ht)
                @info "Success!"
                correct = true
            end
        end
    else
        while !correct
            batchsize = next_batchsize(nprimes_used, batchsize)

            for _ in 1:batchsize
                nprimes_used += 1
                prime = nextluckyprime!(primetracker)
                gens_ff = compute_single_basis(
                    ring, gens_temp_ff, tracer, pairset, ht, 
                    reduced, meta, meta.rng, coeffs_zz, 
                    prime, coeffbuffer
                )
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
        end
    end

    @warn "Used $nprimes_used primes."
    #=
    # TODO
    if meta.usefglm
        fglm_f4!(ring, gens_ff, ht)
    end
    =#

    construct_final_basis(gens_ff, ht, coeffaccum.gb_coeffs_qq)
end

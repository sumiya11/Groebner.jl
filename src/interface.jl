
#######################################

"""
    function groebner(
            polynomials;
            reduced=true,
            ordering=:input,
            certify=false,
            forsolve=false,
            linalg=:exact,
            rng=MersenneTwister(42),
            loglevel=.Warn
    )

Computes a Groebner basis of the ideal generated by array `polynomials`.

If `reduced` is set, returns the reduced basis, which is **unique** (default).

Uses the ordering inherited from `polynomials` by default.
If `ordering` parameter is explicitly specified, it takes precedence.
Possible orderings to specify are

- :input for preserving the input ordering (default)
- :lex for lexicographic
- :deglex for graded lexicographic
- :degrevlex for graded reverse lexicographic

Graded orderings tend to be the fastest.

The algorithm is randomized. Still, by default, the result is **guaranteed to be correct**.

One may set `certify` to `false` to speed up the computation a bit.
The obtained result will be correct with probability > 0.9999

Set `forsolve` to `true` to tell the algorithm to automatically select parameters
for producing a basis that further can used for **solving the input system**. In this case,
the output basis will be in generic position in lexicographic monomial order.
The computation will, however, fail, if the input polynomial system is not solvable.

The `linalg` parameter is responsible for linear algebra backend to be used.
Currently, available options are

- `:exact` for exact sparse linear algebra (default)
- `:prob` for probabilistic sparse linear algebra. Tends to be faster

The algorithm randomized hashing and linear algebra backend depend on random generator `rng`.

Verboseness can be tweaked with the `loglevel` parameter (default is that only warnings are produced).

"""
function groebner(
            polynomials::Vector{Poly};
            reduced::Bool=true,
            ordering::Symbol=:input,
            certify::Bool=false,
            forsolve::Bool=false,
            linalg::Symbol=:exact,
            rng::Rng=Random.MersenneTwister(42),
            loglevel::Logging.LogLevel=Logging.Warn
            ) where {Poly, Rng<:Random.AbstractRNG}

    #= set the logger =#
    prev_logger = Logging.global_logger(ConsoleLogger(stderr, loglevel))

    #= extract ring information, exponents and coefficients
       from input polynomials =#
    # Copies input, so that polynomials would not be changed itself.
    ring, exps, coeffs = convert_to_internal(polynomials, ordering)

    #= check and set algorithm parameters =#
    metainfo = set_metaparameters(ring, ordering, certify, forsolve, linalg)
    # now ring stores computation ordering
    # metainfo is now a struct to store target ordering

    #= change input ordering if needed =#
    assure_ordering!(ring, exps, coeffs, metainfo)

    #= compute the groebner basis =#
    if ring.ch != 0
        # bexps, bcoeffs = groebner_ff(ring, exps, coeffs, reduced, rng, metainfo)
        # if finite field
        # Always returns UInt coefficients #
        bexps, bcoeffs = groebner_ff(ring, exps, coeffs, reduced, rng, metainfo)
    else
        # if rational coefficients
        # Always returns rational coefficients #
        bexps, bcoeffs = groebner_qq(ring, exps, coeffs, reduced, rng, metainfo)
    end

    # ordering in bexps here matches target ordering in metainfo

    #= revert logger =#
    Logging.global_logger(prev_logger)

    # ring contains ordering of computation, it is the requested ordering
    #= convert result back to representation of input =#
    convert_to_output(ring, polynomials, bexps, bcoeffs)
end

#######################################
# Generic isgroebner


"""
    function isgroebner(
                polynomials;
                ordering=:input,
                certify=true,
                rng=MersenneTwister(42),
                loglevel=Logging.Warn
    )

Checks if `polynomials` forms a Groebner basis.

Uses the ordering on `polynomials` by default.
If `ordering` is explicitly specialized, it takes precedence.

The algorithm is randomized. Still, by default, the result is **guaranteed to be correct**.

One may set `certify` to `false` to speed up the computation a bit.
The obtained result will be correct with probability > 0.9999

"""
function isgroebner(
            polynomials::Vector{Poly};
            ordering=:input,
            certify::Bool=true,
            rng::Rng=Random.MersenneTwister(42),
            loglevel::LogLevel=Logging.Warn
    ) where {Poly, Rng<:Random.AbstractRNG}

    #= set the logger =#
    prev_logger = Logging.global_logger(ConsoleLogger(stderr, loglevel))

    #= extract ring information, exponents and coefficients
       from input polynomials =#
    # Copies input, so that polys would not be changed itself.
    ring, exps, coeffs = convert_to_internal(polynomials, ordering)

    #= check and set algorithm parameters =#
    metainfo = set_metaparameters(ring, ordering, certify, false, :exact)
    # now ring stores computation ordering
    # metainfo is now a struct to store target ordering

    #= change input ordering if needed =#
    assure_ordering!(ring, exps, coeffs, metainfo)

    #= compute the groebner basis =#
    if ring.ch != 0
        # if finite field
        # Always returns UInt coefficients #
        flag = isgroebner_ff(ring, exps, coeffs, rng, metainfo)
    else
        # if rational coefficients
        # Always returns rational coefficients #
        flag = isgroebner_qq(ring, exps, coeffs, rng, metainfo)
    end

    #=
    Assuming ordering of `bexps` here matches `ring.ord`
    =#

    #= revert logger =#
    Logging.global_logger(prev_logger)

    flag
end

#######################################
# Generic normalform

"""
    function normalform(
                basis::Vector{Poly},
                tobereduced::Poly;
                ordering::Symbol=:input,
                rng::Rng=Random.MersenneTwister(42),
                loglevel=Logging.Warn
    )

Computes the normal form of `tobereduced` w.r.t `basis`.
The latter is assumed to be a Groebner basis.

Uses the ordering on `basispolys` by default.
If `ordering` is explicitly specialized, it takes precedence.

"""
function normalform(
            basispolys::Vector{Poly},
            tobereduced::Poly;
            ordering::Symbol=:input,
            rng::Rng=Random.MersenneTwister(42),
            loglevel::LogLevel=Logging.Warn
            ) where {Poly, Rng<:Random.AbstractRNG}

    first(normalform(
            basispolys, [tobereduced],
            ordering=ordering, rng=rng, loglevel=loglevel)
    )
end

function normalform(
            basispolys::Vector{Poly},
            tobereduced::Vector{Poly};
            ordering::Symbol=:input,
            rng::Rng=Random.MersenneTwister(42),
            loglevel::LogLevel=Logging.Warn
            ) where {Poly, Rng<:Random.AbstractRNG}

    #= set the logger =#
    prev_logger = Logging.global_logger(ConsoleLogger(stderr, loglevel))

    #= extract ring information, exponents and coefficients
       from input basis polynomials =#
    # Copies input, so that polys would not be changed itself.
    ring1, basisexps, basiscoeffs = convert_to_internal(basispolys, ordering)
    ring2, tbrexps, tbrcoeffs = convert_to_internal(tobereduced, ordering)

    @assert ring1.nvars == ring2.nvars && ring1.ch == ring2.ch
    @assert ring1.ord == ring2.ord

    ring = ring1

    #= check and set algorithm parameters =#
    metainfo = set_metaparameters(ring, ordering, false, false, :exact)

    #= change input ordering if needed =#
    assure_ordering!(ring, basisexps, basiscoeffs, metainfo)
    assure_ordering!(ring, tbrexps, tbrcoeffs, metainfo)


    # We assume basispolys is already a Groebner basis! #

    #= compute the groebner basis =#
    bexps, bcoeffs = normal_form_f4(
                        ring, basisexps, basiscoeffs,
                        tbrexps, tbrcoeffs, rng)

    #=
    Assuming ordering of `bexps` here matches `ring.ord`
    =#

    #= revert logger =#
    Logging.global_logger(prev_logger)

    # ring contains ordering of computation, it is the requested ordering
    #= convert result back to representation of input =#
    convert_to_output(ring, tobereduced, bexps, bcoeffs)
end

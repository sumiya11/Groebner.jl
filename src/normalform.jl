

"""
    function normalform(
                basispolys::Vector{Poly},
                tobereduced::Poly;
                ordering::Symbol=:input,
                randomized::Bool=true,
                rng::Rng=Random.MersenneTwister(42),
                loglevel::Logging.LogLevel=Logging.Warn
    ) where {Poly, Rng<:Random.AbstractRNG}

Computes the normal form of `tobereduced` w.r.t `basispolys`.
The latter is assumed to be a Groebner basis.

Uses the ordering on `basispolys` by default.
If `ordering` is explicitly specialized, it takes precedence.
**(On the fly ordering change not implemented yet :D)**

The algorithm is randomized by default, but
this can be changed with the `randomized` param.
**(Derandomized version not implemented yet :D)**

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

    @assert ring1.nvars == ring2.nvars && ring1.ch == ring2.ch && ring1.ord == ring2.ord
    ring = ring1
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

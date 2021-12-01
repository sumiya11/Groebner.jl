
# The file contains the main algorithm used to compute Groebner bases.
# It composes modular algorithms, fglm, and f4


"""
    Computes the Groebner basis of the ideal generated by the given polynoms

    Parameters:
        . zerodim - if the input ideal is zero-dimensional
                  (i.e. originates from a system with finitely many solutions)
        . reduced - if The reduced Groebner basis to be returned

    Supported input orderings:
        . lex
        . degrevlex

    Invariants:
        . Returns in the same ring the input is given
        .
"""
function groebner(
            fs::Vector{MPoly{Rational{BigInt}}};
            zerodim=false,
            reduced=true)

    Qring = parent(first(fs))

    # determine orderings
    initial_ordering = ordering(Qring)
    @assert initial_ordering in (:lex, :degrevlex)
    if initial_ordering == :lex && zerodim || initial_ordering == :degrevlex
        computing_ordering = :degrevlex
    else
        @warn "Fallback to computations in :lex. This could be slow."
        computing_ordering = :lex
    end
    usefglm = initial_ordering == :lex && computing_ordering == :degrevlex

    # convert polynomials into computing ordering
    fs_ord = change_ordering(fs, computing_ordering)

    println(fs_ord, computing_ordering)
    # scale polynomials to integer coefficients
    fs_ord_zz = scale_denominators(fs_ord)
    Zring = parent(first(fs_ord_zz))

    # reduction moduli
    moduli = Int[ 2^30 + 3 ]
    gbs_gf = [ ]

    i = 0
    while true
        prime = last(moduli)

        # compute the image of fs in GF(prime),
        # by coercing each coefficient into the finite field
        fs_ord_gf = reduce_modulo(fs_ord_zz, prime)

        # groebner basis of the reduced ideal
        gb_ord_gf = f4(fs_ord_gf, reduced=reduced)

        # return back to the initial ordering
        if usefglm
            gb_ord_gf = fglm(gb_ord_gf)
        end

        push!(gbs_gf, gb_ord_gf)

        # TODO: add majority rule based choice

        # trying to reconstruct gbs into rationals
        gb_zz_crt, modulo = reconstruct_crt(gbs_gf, moduli, Zring)
        gb_qq = reconstruct_modulo(gb_zz_crt, modulo, Qring)

        if correctness_checks(gb_qq)
            return reducegb(gb_qq)
        end

        push!(moduli, nextprime(prime))

        i += 1
        if i > 10
            @error "Something probably went wrong in reconstructing in groebner.."
            return
        end
    end

end

function groebner(fs::Vector{MPoly{GFElem{Int}}})
    f4(fs)
end

function correctness_checks(gb_reconstructed)
    return true
end

#------------------------------------------------------------------------------

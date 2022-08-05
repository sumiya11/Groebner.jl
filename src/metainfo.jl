
#------------------------------------------------------------------------------

#=
    GroebnerMetainfo
    Stores parameters for groebner algorithms

    

=#
struct GroebnerMetainfo{Rng}
    # if set, then use fglm algorithm for order conversion
    usefglm::Bool
    # output polynomials monomial order
    targetord::Symbol
    # monomial order for computation
    computeord::Symbol

    # correctness checks levels
    heuristiccheck::Bool
    randomizedcheck::Bool
    guaranteedcheck::Bool

    # linear algebra backend to be used
    # Currently available are
    #   :exact for exact linear algebra
    #   :prob for probabilistic linear algebra
    linalg::Symbol
    
    ground::Symbol

    rng::Rng
end

function set_metaparameters(ring, ordering, certify, forsolve, linalg, rng)
    usefglm = false
    targetord = :lex
    computeord = :lex

    if forsolve
        targetord = :lex
        usefglm = true
        if ordering in (:deglex, :degrevlex)
            computeord = ordering
        else
            computeord = :degrevlex
        end
        # TODO: lex fglm
        computeord = :lex
    else
        if ordering == :input
            ordering = ring.ord
        end
        targetord = ordering
        usefglm = false
        computeord = targetord
    end

    heuristiccheck = true
    # heuristiccheck = false
    randomizedcheck = true
    if certify
        guaranteedcheck = true
    else
        guaranteedcheck = false
    end

    ground = :qq
    if ring.ch == 0
        ground = :qq
    else
        ground = :ff
    end

    @info "Computing in $computeord order, result is in $targetord order"
    @info "Using fglm: $usefglm"

    GroebnerMetainfo(usefglm, targetord, computeord,
                        heuristiccheck, randomizedcheck, guaranteedcheck,
                        linalg, ground, rng)
end

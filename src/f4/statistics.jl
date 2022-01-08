
#=
    The file for keeping track of intermediate results.
    It mirrors some f4 functions to collect statistics
    while computing the basis.
=#

#------------------------------------------------------------------------------

# Contains information about one iteration
struct IterInfo
    basissize::Int
    pairsetsize::Int
    pairsselected::Int
    matrixsize::Tuple{Int, Int}
    matrixdensity::Float64
    newpivots::Int
end

# Stores all collected statistics for one groebner computation
mutable struct Statistics

    data::Vector{IterInfo}

    function Statistics()
        new(Vector{IterInfo}())
    end
end

#------------------------------------------------------------------------------

function f4()

end

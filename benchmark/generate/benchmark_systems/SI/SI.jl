import AbstractAlgebra

include((@__DIR__) * "/for_gleb.jl")
include((@__DIR__) * "/parser.jl")

function load_SI_problem(name)
    if name == "SIWR"
        sys = read_SIWR()[2]
        @assert AbstractAlgebra.ordering(AbstractAlgebra.parent(sys[1])) == :degrevlex
        return sys
    elseif name == "SEAIJRC"
        sys = read_SEAIJRC()[2]
        @assert AbstractAlgebra.ordering(AbstractAlgebra.parent(sys[1])) == :degrevlex
        return sys
    else
        throw("Beda beda ogorchenie")
    end
end

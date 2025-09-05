import AbstractAlgebra

include((@__DIR__) * "/for_gleb.jl")
include((@__DIR__) * "/parser.jl")

function load_SI_problem(name; k=AbstractAlgebra.QQ)
    if name == "SIWR"
        sys = read_SIWR(k=k)[2]
        @assert AbstractAlgebra.internal_ordering(AbstractAlgebra.parent(sys[1])) == :degrevlex
        return sys
    elseif name == "SEAIJRC"
        sys = read_SEAIJRC(k=k)[2]
        @assert AbstractAlgebra.internal_ordering(AbstractAlgebra.parent(sys[1])) == :degrevlex
        return sys
    else
        throw("Beda beda ogorchenie")
    end
end

import AbstractAlgebra

include((@__DIR__) * "/parser.jl")

function load_MQ_problem(name)
    sys = read_MQ_GF(name)
    @assert AbstractAlgebra.internal_ordering(AbstractAlgebra.parent(sys[1])) == :degrevlex
    sys
end

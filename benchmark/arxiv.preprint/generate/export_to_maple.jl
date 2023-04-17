using Groebner
using AbstractAlgebra

include((@__DIR__)*"/benchmark_systems.jl")

header = """
with(Groebner):
with(PolynomialIdeals):

"""

template = """
J := PolynomialIdeal({(1)}, charactesistic=(3)):
print("Running (4)");
st := time[real]():
Groebner[Basis](J, tdeg((2)), method=direct):
print("(4): ", time[real]() - st);
"""

function export_maple(io, name, system)
    global ground
    R = parent(first(system))
    t1 = join(map(string, system), ", ")
    t2 = join(map(string, gens(R)), ", ")
    t3 = string(characteristic(base_ring(R)))
    t4 = name
    s = replace(template, "(1)"=>t1,"(2)"=>t2,"(3)"=>t3,"(4)"=>t4)
    println(io, s)
end

function generate(flag)
    if flag
        ground = GF(2^31-1)
        systems = benchmark_systems_ff(ground)
    else
        ground = QQ
        systems = benchmark_systems_qq(ground)
    end
    
    io = open((@__DIR__)*"/generated_maple_$(flag ? "ff" : "qq").mpl", "w")
    
    println(io, header)

    for (name, system) in systems
        export_maple(io, name, system)
    end
    
    close(io) 
end

generate(true)
generate(false)

using Groebner
using AbstractAlgebra

header = """
with(Groebner):
with(PolynomialIdeals):

kernelopts(numcpus=1);

"""

template = """
J := [(1)]:
print("Running (4)");
st := time[real]():
Groebner[Basis](J, tdeg((2)), method=fgb, characteristic=(3)):
print("(4): ", time[real]() - st);
"""

function export_maple(io, name, system)
    global ground
    R = parent(first(system))
    t1 = join(map(string, system), ", ")
    t2 = join(map(string, gens(R)), ", ")
    t3 = string(characteristic(base_ring(R)))
    t4 = name
    s = replace(template, "(1)" => t1, "(2)" => t2, "(3)" => t3, "(4)" => t4)
    println(io, s)
end

function generate()
    ground = GF(2^31 - 1)
    systems = [
        ("cylic7", Groebner.cyclicn(7, ground=ground)),
        ("cylic8", Groebner.cyclicn(8, ground=ground)),
        ("katsura10", Groebner.katsuran(10, ground=ground)),
        ("katsura11", Groebner.katsuran(11, ground=ground)),
        ("eco12", Groebner.eco12(ground=ground)),
        ("eco13", Groebner.eco13(ground=ground)),
        ("noon7", Groebner.noonn(7, ground=ground)),
        ("noon8", Groebner.noonn(8, ground=ground))
    ]

    io = open((@__DIR__) * "/generated_maple.mpl", "w")

    println(io, header)

    for (name, system) in systems
        export_maple(io, name, system)
    end

    close(io)
end

generate()

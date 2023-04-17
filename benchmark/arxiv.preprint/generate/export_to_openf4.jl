
using AbstractAlgebra
using Groebner

include((@__DIR__)*"/benchmark_systems.jl")

template = "
int (5)F4(bool magma)
{
    cout << \"#########################################################\" << endl;
    cout << \"#                         (1)                        #\" << endl;
    cout << \"#########################################################\" << endl << endl;
    eltType::setModulo(modulo);
    int nbGen;
    Monomial::initMonomial((2));
    vector<Polynomial<eltType>> (3);
    (4)
    Ideal<eltType> (5)((3), (2), 1000000);
    nbGen=(5).f4();
    if(magma)
    {
        (5).printReducedGroebnerBasis(\"(1)\", modulo);
    }
    return nbGen;
}
"

function export_openf4(io, name, system)
    global ground
    nv = nvars(parent(first(system)))
    newR, _ = PolynomialRing(ground, ["x$i" for i in 0:nv-1], ordering=:degrevlex)
    system = map(f -> change_base_ring(ground, f, parent=newR), system)
    t1 = uppercase(name)
    t2 = string(nv)
    t3 = "pol$name"
    t4 = "\n\t////////////\n"
    t4 *= "\t//\t\t$(uppercase(name))\t\t//\n"
    t4 *= "\t// $(nvars(parent(first(system)))) vars: $(gens(parent(first(system))))\n"
    for poly in system
        t4 *= "\tpol$name.emplace_back(\"$poly\");\n"
    end
    t4 *= "\t////////////\n\n"
    t5 = name
    s = replace(template, "(1)"=>t1,"(2)"=>t2,"(3)"=>t3,"(4)"=>t4,"(5)"=>t5)
    println(io, s)
end

ground = GF(2^31-1)
systems = benchmark_systems_ff(ground)

io = open((@__DIR__)*"/../openf4_systems.txt", "w")

for (name, system) in systems
    export_openf4(io, name, system)
end

close(io)

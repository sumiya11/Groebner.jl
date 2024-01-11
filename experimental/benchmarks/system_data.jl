
using AbstractAlgebra
using Base.Threads

# import Singular

include("../src/Groebner.jl")

function generatesystems(ground)
    [
        # ("noon 3", Groebner.noon3(k=ground)),
        ("cyclic 4", Groebner.cyclicn(4, k=ground)),
        ("cyclic 5", Groebner.cyclicn(5, k=ground)),
        ("cyclic 6", Groebner.cyclicn(6, k=ground)),
        ("cyclic 7", Groebner.cyclicn(7, k=ground)),
        ("cyclic 8", Groebner.cyclicn(8, k=ground)),
        ("cyclic 9", Groebner.cyclicn(9, k=ground))
        # ("eco 11", Groebner.eco11(k=ground)),
        # ("eco 12", Groebner.eco12(k=ground)),
        # ("eco 13", Groebner.eco13(k=ground))
        #("noon 6", Groebner.noonn(6, k=ground)),
        #("noon 7", Groebner.noonn(7, k=ground)),
        #("noon 8", Groebner.noonn(8, k=ground)),
        # ("noon 9", Groebner.noonn(9, k=ground)),
        #("cyclic 10", Groebner.rootn(10, k=ground)),
        #("cyclic 11", Groebner.rootn(11, k=ground)),
        #("cyclic 12", Groebner.rootn(12, k=ground)),
        #("henrion 5", Groebner.henrion5(k=ground)),
        #("henrion 6", Groebner.henrion6(k=ground)),
        # ("henrion 7", Groebner.henrion7(k=ground))
    ]
end

function iszerodim(gb)
    R = parent(first(gb))
    vs = gens(R)
    for v in vs
        flag = false
        for lt in map(leading_monomial, gb)
            if exponent_vector(v, 1) .* total_degree(lt) == exponent_vector(lt, 1)
                flag = true
            end
        end
        if !flag
            return false
        end
    end
    return true
end

function nsols(gb)
    R = parent(first(gb))
    vs = gens(R)
    sols = 1
    for v in vs
        for lt in map(leading_monomial, gb)
            if exponent_vector(v, 1) .* total_degree(lt) == exponent_vector(lt, 1)
                sols *= total_degree(lt)
            end
        end
    end
    return sols
end

#=
function systeminfo(system)
    R = AbstractAlgebra.parent(system[1])
    modulo = AbstractAlgebra.characteristic(R)
    n = AbstractAlgebra.nvars(R)
    ground_s = Singular.N_ZpField(modulo)
    R_s, _ = Singular.polynomial_ring(ground_s, ["x$i" for i in 1:n], ordering=:degrevlex)

    system_s = map(
        f -> AbstractAlgebra.change_base_ring(
                    ground_s,
                    AbstractAlgebra.map_coefficients(c -> ground_s(c.d), f),
                    parent=R_s),
        system)

    ideal_s = Singular.Ideal(R_s, system_s)
    gb = Singular.std(ideal_s)

    # println("Is zerodim ", Singular.is_zerodim(ideal_s))
    println("dim R / I ", length(Singular.gens(Singular.kbase(gb))))
end
=#

function systeminfo(system)
    println(iszerodim(Groebner.groebner(system, ordering=:degrevlex)))
end

function runall(ground)
    println("-"^20)
    for (name, system) in generatesystems(ground)
        println(name)
        systeminfo(system)
        println("-"^20)
    end
end

ground = GF(2^31 - 1)
runall(ground)

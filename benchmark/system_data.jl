
using AbstractAlgebra
using Base.Threads


function generatesystems(ground)
    [
        ("noon 3", Groebner.noon3(ground=ground)),
        # ("katsura 11", Groebner.katsuran(11, ground=ground)),
        # ("katsura 12", Groebner.katsuran(12, ground=ground)),
        # ("katsura 13", Groebner.katsuran(13, ground=ground))
        # ("eco 11", Groebner.eco11(ground=ground)),
        # ("eco 12", Groebner.eco12(ground=ground)),
        # ("eco 13", Groebner.eco13(ground=ground))
        # ("noon 7", Groebner.noonn(7, ground=ground)),
        # ("noon 8", Groebner.noonn(8, ground=ground)),
        # ("noon 9", Groebner.noonn(9, ground=ground)),
        # ("cyclic 12", Groebner.rootn(12, ground=ground)),
        # ("cyclic 13", Groebner.rootn(13, ground=ground)),
        # ("cyclic 14", Groebner.rootn(14, ground=ground)),
        ("henrion 5", Groebner.henrion5(ground=ground)),
        ("henrion 6", Groebner.henrion6(ground=ground)),
        ("henrion 7", Groebner.henrion7(ground=ground))
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

function systemdata(ground)

    println("-"^20)
    for (name, system) in generatesystems(ground)
        println(name)
        gb = Groebner.groebner(system, ordering=:degrevlex)
        # println(gb)
        println("# LEAD = ", length(leading_monomial.(gb)))
        println("ZERODIM? ", iszerodim(gb))
        if iszerodim(gb)
            println("NSOLS ", nsols(gb))
        end
        println("-"^20)
    end
end

ground = GF(2^31-1)
# systemdata(ground)

function uwu(x)
    sleep(rand())
    println("from $(threadid()) print $x")
end

function runall()
    for i in 1:10
        task = @spawn uwu(i)
        wait(task)
    end
end


runall()

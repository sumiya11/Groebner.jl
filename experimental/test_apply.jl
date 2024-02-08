using AbstractAlgebra
using Primes
using Random
include("../src/Groebner.jl")

function fooo()
    R, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering=:degrevlex)

    # sys = [x * y + y, x * y + x + y]
    # sys = [x * y^3 + y, x * y + x + y]
    sys = Groebner.katsuran(9)
    R = parent(sys[1])
    println(sys)
    Groebner.groebner(sys)

    Random.seed!(42)
    A, B = 0, 0

    boot = 100
    for p in Primes.nextprimes(2^25, 100)
        @info "" p boot
        for i in 1:boot
            if (i % 10) == 0
                @info "$i / $boot"
            end
            #
            # sys_x = map(f -> evaluate(f, xy .* gens(R)), sys)
            sys_x = empty(sys)
            for f in sys
                _f = zero(f)
                for t in monomials(f)
                    __c = rand(-100:100)
                    if __c == 0
                        __c = 1
                    end
                    _f += t * __c
                end
                push!(sys_x, _f)
            end
            #
            sys_x_mod_p = map(f -> map_coefficients(c -> GF(p)(numerator(c)), f), sys_x)
            if all(iszero, sys_x_mod_p)
                continue
            end
            trace, gb = Groebner.groebner_learn(sys_x_mod_p)
            backup = deepcopy(trace)
            for p2 in Primes.nextprimes(2^25, 10)
                if p == p2
                    continue
                end
                # trace, gb = Groebner.groebner_learn(sys_x_mod_p)
                sys_x_mod_p2 =
                    map(f -> map_coefficients(c -> GF(p2)(numerator(c)), f), sys_x)
                # @info "" p p2
                # println(sys_x)
                # println(sys_x_mod_p)
                # println(sys_x_mod_p2)
                # println(Groebner.groebner(sys_x_mod_p))
                # println(Groebner.groebner(sys_x_mod_p2))
                success, gb2 = Groebner.groebner_apply!(trace, sys_x_mod_p2)
                if !success
                    println(A / B)
                    A += 1
                    # @info "$success $p $p2"
                    trace = backup
                    backup = deepcopy(trace)
                end
                B += 1
            end
        end
    end
    A / B
end

fooo()

#=
p,p2 = 3,19

sys =           [-11*x*y + 53*y, 83*x*y + x - 70*y]
sys_mod_p =     [x*y + 2*y, 2*x*y + x + 2*y]
sys_mod_p2 =    [8*x*y + 15*y, 7*x*y + x + 6*y]
gb_mod_p =      [x + y, y^2 + y]
gb_mod_p2 =     [y, x]
=#

R, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering=:degrevlex)
p, p2 = 3, 19

sys = [-11 * x * y + 53 * y, 83 * x * y + x - 70 * y]
sys_mod_p = map(
    f -> map_coefficients(c -> GF(p)(numerator(c)), f),
    [x * y + 2 * y, 2 * x * y + x + 2 * y]
)
sys_mod_p2 = map(
    f -> map_coefficients(c -> GF(p2)(numerator(c)), f),
    [8 * x * y + 15 * y, 7 * x * y + x + 6 * y]
)

Groebner.groebner(sys_mod_p, loglevel=-100)
Groebner.groebner(sys_mod_p2, loglevel=-100)

trace, gb = Groebner.groebner_learn(sys_mod_p)
flag, gb = Groebner.groebner_apply!(trace, sys_mod_p, loglevel=-5)
flag, gb = Groebner.groebner_apply!(trace, sys_mod_p2, loglevel=-5)

trace, gb = Groebner.groebner_learn(sys_mod_p2)
flag, gb = Groebner.groebner_apply!(trace, sys_mod_p2)
flag, gb = Groebner.groebner_apply!(trace, sys_mod_p)

###############

p, p2 = nextprime(2^27), nextprime(2^27)

sys = Groebner.cyclicn(8)

sys_mod_p = map(f -> map_coefficients(c -> GF(p)(numerator(c)), f), sys)
sys_mod_p2 = map(f -> map_coefficients(c -> GF(p2)(numerator(c)), f), sys)

a = length.(Groebner.groebner(sys_mod_p))
b = length.(Groebner.groebner(sys_mod_p2))
a == b

trace, gb = Groebner.groebner_learn(sys_mod_p);
flag, gb = Groebner.groebner_apply!(trace, sys_mod_p);
flag, gb = Groebner.groebner_apply!(trace, sys_mod_p2);

@profview Groebner.groebner_apply!(trace, sys_mod_p)

trace, gb = Groebner.groebner_learn(sys_mod_p2)
flag, gb = Groebner.groebner_apply!(trace, sys_mod_p2)
flag, gb = Groebner.groebner_apply!(trace, sys_mod_p)

#=
p,p2 = 3,11

sys =           [-x*y^3 + 19*y, 67*x*y + 100*x - 73*y]
sys_mod_p =     [2*x*y^3 + y, x*y + x + 2*y]
sys_mod_p2 =    [10*x*y^3 + 8*y, x*y + x + 4*y]
gb_mod_p =      [x*y + x + 2*y]
gb_mod_p2 =     [x*y + x + 4*y, y^3 + 10*y^2 + 3*x + 3*y]
=#

R, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering=:degrevlex)
p, p2 = 3, 13

sys = [-x * y^3 + 19 * y, 67 * x * y + 100 * x - 73 * y]
sys_mod_p = map(
    f -> map_coefficients(c -> GF(p)(numerator(c)), f),
    [2 * x * y^3 + y, x * y + x + 2 * y]
)
sys_mod_p2 = map(
    f -> map_coefficients(c -> GF(p2)(numerator(c)), f),
    [10 * x * y^3 + 8 * y, x * y + x + 4 * y]
)

Groebner.groebner(sys_mod_p, loglevel=-100)
Groebner.groebner(sys_mod_p2, loglevel=-100)

trace, gb = Groebner.groebner_learn(sys_mod_p)
flag, gb = Groebner.groebner_apply!(trace, sys_mod_p)
flag, gb = Groebner.groebner_apply!(trace, sys_mod_p2)

#=
┌ Info: 
│   p = 3
└   p2 = 31
sys =           [4*x*y^3 - 77*y, 98*x*y + 25*x - 34*y]
sys_mod_p =     [x*y^3 + y, 2*x*y + x + 2*y]
sys_mod_p2 =    [4*x*y^3 + 16*y, 5*x*y + 25*x + 28*y]
gb_mod_p =      [x*y + 2*x + y]
gb_mod_p2 =     [x*y + 5*x + 18*y, x^2 + 8*y^2 + 15*x + 13*y, y^3 + 26*y^2 + 19*x + 11*y]
=#

R, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering=:degrevlex)
p, p2 = 3, 31

sys = [4 * x * y^3 - 77 * y, 98 * x * y + 25 * x - 34 * y]
sys_mod_p = map(
    f -> map_coefficients(c -> GF(p)(numerator(c)), f),
    [x * y^3 + y, 2 * x * y + x + 2 * y]
)
sys_mod_p2 = map(
    f -> map_coefficients(c -> GF(p2)(numerator(c)), f),
    [4 * x * y^3 + 16 * y, 5 * x * y + 25 * x + 28 * y]
)

Groebner.groebner(sys_mod_p, loglevel=-100)
Groebner.groebner(sys_mod_p2, loglevel=-100)

trace, gb = Groebner.groebner_learn(sys_mod_p)
flag, gb = Groebner.groebner_apply!(trace, sys_mod_p)
flag, gb = Groebner.groebner_apply!(trace, sys_mod_p2)

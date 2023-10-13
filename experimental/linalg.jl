using Groebner
using AbstractAlgebra
using PrettyTables
using ProgressMeter
using Primes

fields = vcat(
    [GF(p) for p in Primes.primes(2, 200)],
    [GF(p) for p in Primes.primes(201, 500)[1:3:end]],
    [GF(p) for p in Primes.primes(501, 1000)[1:5:end]],
    [GF(p) for p in Primes.primes(1001, 2000)[1:15:end]],
    [GF(p) for p in Primes.primes(5001, 10000)[1:40:end]]
    # [GF(p) for p in Primes.primes(10_000, 50_000)[1:150:end]],
    # [GF(p) for p in Primes.primes(50_001, 300_000)[1:500:end]],
    # [GF(p) for p in Primes.primes(300_001, 2^20)[1:1000:end]]
)

progress = Progress(
    length(fields),
    dt=0.1,
    barglyphs=BarGlyphs('|', '█', ['▁', '▂', '▃', '▄', '▅', '▆', '▇'], ' ', '|')
    # barlen=10
)
names = [
    :noon2,
    :noon3,
    :noon4,
    :noon5,
    :eco5,
    :cycl2,
    :cycl3,
    :cycl4,
    :cycl5,
    :kat2,
    :kat3,
    :kat4,
    :kat5,
    :kat6
]

m, n = length(names), length(fields)
boot = 10
table = Array{Any, 2}(undef, n, m)
for (i, field) in enumerate(fields)
    R, (x, y) = PolynomialRing(field, ["x", "y"])
    cases = [
        (Groebner.noonn(2, ground=field, ordering=:degrevlex)),
        (Groebner.noonn(3, ground=field, ordering=:degrevlex)),
        (Groebner.noonn(4, ground=field, ordering=:degrevlex)),
        (Groebner.noonn(5, ground=field, ordering=:degrevlex)),
        (Groebner.eco5(ground=field, ordering=:degrevlex)),
        (Groebner.cyclicn(2, ground=field, ordering=:degrevlex)),
        (Groebner.cyclicn(3, ground=field, ordering=:degrevlex)),
        (Groebner.cyclicn(4, ground=field, ordering=:degrevlex)),
        (Groebner.cyclicn(5, ground=field, ordering=:degrevlex)),
        (Groebner.katsuran(2, ground=field, ordering=:degrevlex)),
        (Groebner.katsuran(3, ground=field, ordering=:degrevlex)),
        (Groebner.katsuran(4, ground=field, ordering=:degrevlex)),
        (Groebner.katsuran(5, ground=field, ordering=:degrevlex)),
        (Groebner.katsuran(6, ground=field, ordering=:degrevlex))
    ]
    for (j, system) in enumerate(cases)
        flag = true
        xs = gens(parent(system[1]))
        for _ in 1:boot
            gb = groebner(system)
            flag &= isgroebner(gb)
            flag &= all(iszero, normalform(gb, system))
            dilation = map(_ -> rand(field), 1:length(xs))
            shift = map(_ -> rand(field), 1:length(xs))
            system = map(poly -> evaluate(poly, xs .* dilation .+ shift), system)
        end
        table[i, j] = flag
    end
    next!(progress)
end
finish!(progress)

highlighter = Highlighter((v, i, j) -> !v[i, j], crayon"red bold")

pretty_table(
    table,
    limit_printing=false,
    header=names,
    display_size=(-1, -1),
    title="Is basis correct?",
    highlighters=(highlighter,),
    row_labels=map(f -> Symbol(:GF_, Symbol(string(characteristic(f)))), fields)
)
print("Out of $(m*n) tests: ")
failed = m * n - sum(table)
if failed > 0
    printstyled("$failed failed\n", color=:light_red)
else
    printstyled("All passed\n", color=:light_green)
end

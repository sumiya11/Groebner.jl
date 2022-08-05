
import AbstractAlgebra

function read_MQ_GF(filename)
    @info "Reading MQ $filename"
    polys = []
    apath = abspath("") * "benchmark\\MQ\\$filename"
    open(apath, "r") do fio
        line1, line2, line3 = readline(fio), readline(fio), readline(fio)
        ch = parse(Int, line1[findfirst('(', line1)+1:findfirst(')', line1)-1])
        nv = parse(Int, line2[findfirst(':', line2)+1:end])
        ne = parse(Int, line3[findfirst(':', line3)+1:end])

        R, xs = AbstractAlgebra.PolynomialRing(AbstractAlgebra.GF(ch), nv, ordering=:degrevlex)
        @info "Created ring" R xs

        labels = sort!(union([x*y for x in xs for y in xs], xs, [R(1)]), rev=true)

        rest = readlines(fio)
        rest = rest[5:end]

        for (i, s) in enumerate(rest)
            @info "Reading poly $i/$(length(rest))"
            s = replace(s, " ;" => "")
            poly = sum(map(x -> parse(Int, x), split(s, " ")) .* labels)
            push!(polys, poly)
        end
    end
    @info "Read" length(polys)
    R = parent(first(polys))
    map(R, polys)
end

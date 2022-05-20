
import Nemo

function read_SEAIJRC()
    @info "Reading SEAIJRC"
    polys = []
    apath = "C:\\data\\projects\\gbdata\\gbSEAIJRC"
    open(apath, "r") do f
        head = "b, alpha, g2, Ninv, k, g1, q, r, sat_aux1, sat_aux2"
        body = readline(f)
        @info "First line" head
        headv = map(x -> Symbol(x), split(head, ", "))
        e = Meta.parse(join(headv, ","))
        R, xs = eval(:((R, $e) = Nemo.PolynomialRing(Nemo.QQ, $headv)))
        @info "Created ring" R xs
        rrr = split(body, ",")
        for (i, s) in enumerate(rrr)
            @info "Reading poly $i/$(length(rrr))"
            poly = eval(Meta.parse(s))
            push!(polys, poly)
        end
    end
    @info "Read" length(polys)
    R.(polys)
end

function read_SIWR()
    @info "Reading SIWR"
    polys = []
    apath = "C:\\data\\projects\\gbdata\\gbSIWR"
    open(apath, "r") do f
        head = readline(f)
        body = readline(f)
        headv = map(x -> Symbol(x), split(head, ", "))
        e = Meta.parse(join(headv, ","))
        R, xs = eval(:((R, $e) = Nemo.PolynomialRing(Nemo.QQ, $headv)))
        @info "Created ring" R xs
        rrr = split(body, ",")
        for (i, s) in enumerate(rrr)
            @info "Reading poly $i/$(length(rrr))"
            poly = eval(Meta.parse(s))
            push!(polys, poly)
        end
    end
    @info "Read" length(polys)
    R.(polys)
end

function read_MAPK()
    @info "Reading MAPK"
    polys = []
    apath = "../../MAPK"
    # apath = "MAPK"
    open(apath, "r") do f
        head = "c0001, a10, gamma1000, alpha10, b00, beta11, c0111, beta10, alpha11, beta01, alpha01, gamma1100, c0011, c0010, b10, gamma1101, a00, b01, c1011, a01, gamma1110, gamma0100, sat_aux1, sat_aux2, sat_aux3, sat_aux4, sat_aux5"
        body = readline(f)
        @info "First line" head
        headv = map(x -> Symbol(x), split(head, ", "))
        e = Meta.parse(join(headv, ","))
        R, xs = eval(:((R, $e) = Nemo.PolynomialRing(Nemo.QQ, $headv)))
        @info "Created ring" R xs
        rrr = split(body, ",")
        l = length(rrr)
        for i in 1:length(rrr)
            @info "Reading poly $i/$(length(rrr))"
            str = map(strip, split(rrr[1], "+"))
            poly = R(0)
            for s1 in str
                # println(s1)
                # println("_________________")
                if occursin("-", s1)
                    sp1 = split(s1, "-")
                    if startswith(s1, "-")
                        sp1 = sp1[2:end]
                        poly -= eval(Meta.parse(sp1[1]))
                    else
                        poly += eval(Meta.parse(sp1[1]))
                    end
                    for s2 in sp1[2:end]
                        poly -= eval(Meta.parse(s2))
                    end
                else
                    poly += eval(Meta.parse(s1))
                end
            end
            push!(polys, poly)
            popfirst!(rrr)
        end
    end
    @info "Read" length(polys)
    R.(polys)
end

function read_MAPK_include()
    @info "Reading MAPK"
    polys = []
    apath = "../../MAPK"

    head = "c0001, a10, gamma1000, alpha10, b00, beta11, c0111, beta10, alpha11, beta01, alpha01, gamma1100, c0011, c0010, b10, gamma1101, a00, b01, c1011, a01, gamma1110, gamma0100, sat_aux1, sat_aux2, sat_aux3, sat_aux4, sat_aux5"
    @info "First line" head
    headv = map(x -> Symbol(x), split(head, ", "))
    e = Meta.parse(join(headv, ","))
    R, xs = eval(:((R, $e) = Nemo.PolynomialRing(Nemo.QQ, $headv)))

    collect(include(apath))
end

function parse_include(filename)
    head = "x, y, z"
    headv = map(x -> Symbol(x), split(head, ", "))
    e = Meta.parse("x, y, z")
    R, xs = eval(:((R, $e) = Nemo.PolynomialRing(Nemo.QQ, $headv)))

    include(abspath("benchmark\\data\\$filename"))
end

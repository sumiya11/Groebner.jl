import Nemo
using SHA

function extract_ring(line1, line2)
    vars_str = map(f -> string(strip(f, [',', ' ', '\t'])), split(line1, ","))
    char = parse(BigInt, strip(line2))
    @assert char < typemax(UInt)
    base_field = iszero(char) ? Nemo.QQ : Nemo.GF(UInt(char))
    ring_nemo, vars_nemo = Nemo.PolynomialRing(base_field, vars_str, ordering=:degrevlex)
    base_field, ring_nemo, vars_nemo
end

function parse_polynomial_from_terms(
    ring_nemo,
    constant_type::T,
    terms_exploded_str,
    str_to_var_idx
) where {T}
    base_field = Nemo.base_ring(ring_nemo)
    n = Nemo.nvars(ring_nemo)
    cfs = Vector{Nemo.elem_type(base_field)}(undef, length(terms_exploded_str))
    exps = Vector{Vector{Int}}(undef, length(terms_exploded_str))
    for (i, t) in enumerate(terms_exploded_str)
        exp_ = zeros(Int, n)
        cf_ = one(base_field)
        for m in t
            varexp = split(m, "^")
            if length(varexp) == 1
                x = varexp[1]
                if haskey(str_to_var_idx, x)
                    exp_[str_to_var_idx[x]] += 1
                else
                    constant = parse(constant_type, x)
                    cf_ *= base_field(constant)
                end
            else
                @assert length(varexp) == 2
                x, y = varexp
                exp_[str_to_var_idx[x]] += parse(Int, y)
            end
        end
        cfs[i] = cf_
        exps[i] = exp_
    end
    cfs, exps
end

function getlines_backend_dependent(result)
    if startswith(result, "#")
        lines = split(result, "\n")
        line2 = strip(split(lines[1])[end])
        line1 = join(split(lines[2])[4:end], " ")
        lines_polys =
            filter(!isempty, map(f -> string(strip(f, [' ', '\n', '\r'])), lines[5:end]))
        @assert startswith(lines_polys[1], "[")
        lines_polys[1] = lines_polys[1][2:end]
        @assert endswith(lines_polys[end], "]:")
        lines_polys[end] = lines_polys[end][1:(end - 2)]
        return vcat(line1, line2, lines_polys)
    end
    if contains(result, ")") && contains(result, "(")
        result = replace(result, "(" => "", ")" => "")
    end
    filter(!isempty, split(result, "\n"))
end

function parse_system_naive(result::String)
    lines = getlines_backend_dependent(result)
    @assert length(lines) >= 2
    base_field, ring_nemo, vars_nemo = extract_ring(lines[1], lines[2])
    if length(lines) == 2
        return ring_nemo, []
    end
    polys_str = map(s -> string(strip(s, [' ', ',', '\t'])), lines[3:end])
    polys_str = filter(!isempty, polys_str)
    str_to_var = Dict{String, Nemo.elem_type(ring_nemo)}(string(v) => v for v in vars_nemo)
    str_to_var_idx = Dict{String, Int}(string(v) => i for (i, v) in enumerate(vars_nemo))
    polys_nemo = Vector{Nemo.elem_type(ring_nemo)}()
    constant_type = base_field == Nemo.QQ ? Rational{BigInt} : Int64
    for poly_str in polys_str
        terms_str_plus = split(poly_str, '+')
        terms_str = empty(terms_str_plus)
        @assert length(terms_str_plus) > 0
        for term_str in terms_str_plus
            term_str_minus = split(term_str, '-')
            if !isempty(strip(term_str_minus[1]))
                push!(terms_str, strip(term_str_minus[1]))
            end
            append!(terms_str, map(s -> "-1*$(strip(s))", term_str_minus[2:end]))
        end
        terms_exploded_str = map(t -> map(strip, split(t, "*")), terms_str)
        cfs, exps = parse_polynomial_from_terms(
            ring_nemo,
            constant_type,
            terms_exploded_str,
            str_to_var_idx
        )
        poly_nemo = ring_nemo(cfs, exps)
        push!(polys_nemo, poly_nemo)
    end
    ring_nemo, polys_nemo
end

# "inline" tests
for parse_system in [parse_system_naive]
    ring, polys = parse_system("""
    x, y,z  
    13
    17
    x - y + 13*z - 4*x - 3*y + 18*z,
    x^3*y + 14,
    """)

    flag1 = Nemo.symbols(ring) == [:x, :y, :z]
    flag2 = Nemo.base_ring(ring) == Nemo.GF(13)
    flag3 = map(string, polys) == ["4", "10*x + 9*y + 5*z", "x^3*y + 1"]
    @assert flag1 && flag2 && flag3 "Parsing routine $parse_system is broken for groebner"

    ring, polys = parse_system("""
    x1,x2  
    0
    5*x1^2 - 10*x2 + 1267650600228229401496703205376// 1267650600228229401496703205375,
    -x2^20*x1,
    -2*x1 - 3*x2,
    -1,
    0
    """)

    flag1 = Nemo.symbols(ring) == [:x1, :x2]
    flag2 = Nemo.base_ring(ring) == Nemo.QQ
    flag3 =
        map(string, polys) == [
            "5*x1^2 - 10*x2 + 1267650600228229401496703205376//1267650600228229401496703205375",
            "-x1*x2^20",
            "-2*x1 - 3*x2",
            "-1",
            "0"
        ]
    @assert flag1 && flag2 && flag3 "Parsing routine $parse_system is broken for groebner"

    ring, polys = parse_system("""
    #Reduced Groebner basis for input in characteristic 1073741827
    #for variable order z1, z2, z3, z4, z5, z6, z7
    #w.r.t. grevlex monomial ordering
    #consisting of 209 elements:
    [1*z1^1+1*z2^1+1*z3^1+1*z4^1+1*z5^1+1*z6^1+1*z7^1,
    -z6^3]:
    """)
    flag1 = Nemo.symbols(ring) == [:z1, :z2, :z3, :z4, :z5, :z6, :z7]
    flag2 = Nemo.base_ring(ring) == Nemo.GF(2^30 + 3)
    flag3 = map(string, polys) == ["z1 + z2 + z3 + z4 + z5 + z6 + z7", "1073741826*z6^3"]
    @assert flag1 && flag2 && flag3 "Parsing routine $parse_system is broken for msolve"

    ring, polys = parse_system("""
    x, y  
    13
    (1*x^1) + (1*y^1),
    (1*y^2) + (16*1),
    """)
    flag1 = Nemo.symbols(ring) == [:x, :y]
    flag2 = Nemo.base_ring(ring) == Nemo.GF(13)
    flag3 = map(string, polys) == ["x + y", "y^2 + 3"]
    @assert flag1 && flag2 && flag3 "Parsing routine $parse_system is broken for openf4"
end

function hash_polynomial(poly, h)
    poly_str = repr(poly)
    sha3_256(poly_str)
end

function get_certificate(ring, polys)
    filter!(!iszero, polys)
    @assert !isempty(polys)
    sort!(polys, by=Nemo.leading_monomial)
    polys = map(poly -> Nemo.divexact(poly, Nemo.leading_coefficient(poly)), polys)
    h = map(i -> UInt8(0), 1:32)
    for (i, poly) in enumerate(polys)
        i_8 = i % UInt8
        h = i_8 .⊻ hash_polynomial(poly, h)
    end
    certificate = bytes2hex(h)
    certificate
end

function is_certificate_standardized(certificate)
    certificate == standardize_certificate(certificate)
end

function standardize_certificate(certificate)
    certificate = strip(string(certificate), ['\n', ' ', '\r'])
    certificate
end

function compute_basis_validation_hash(result)
    ring, polys = parse_system_naive(result)
    if isempty(polys)
        return false, ""
    end
    certificate = get_certificate(ring, polys)
    true, standardize_certificate(certificate)
end

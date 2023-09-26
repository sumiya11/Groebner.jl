import Nemo

function parse_system_using_meta_parse(result::String)
    lines = split(result, "\n")
    @assert length(lines) > 2
    base_field, ring_nemo, vars_nemo = extract_ring(lines[1], lines[2])
    polys_str = map(s -> strip(s, [' ', ',', '\t']), lines[3:end])
    polys_str = filter(!isempty, polys_str)
    polys_nemo = Vector{Nemo.elem_type(ring_nemo)}()
    var_mapping = Dict{Symbol, Nemo.elem_type(ring_nemo)}(
        x => gen(ring_nemo, i) for (i, x) in enumerate(symbols(ring_nemo))
    )
    @info "" polys_str var_mapping ring_nemo
    for poly_str in polys_str
        poly_expr = Meta.parse(poly_str)
        poly_nemo = myeval(poly_expr, base_field, var_mapping)
        push!(polys_nemo, poly_nemo)
    end
    ring_nemo, polys_nemo
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

function parse_system_naive(result::String)
    lines = split(result, "\n")
    @assert length(lines) > 2
    base_field, ring_nemo, vars_nemo = extract_ring(lines[1], lines[2])
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
for parse_system in [parse_system_naive] #parse_system_using_meta_parse]
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
    @assert flag1 && flag2 && flag3 "Parsing routine $parse_system is broken"

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
    @assert flag1 && flag2 && flag3 "Parsing routine $parse_system is broken"
end

function get_certificate(ring, polys)
    length(polys)
end

function is_certificate_standardized(certificate)
    certificate == standardize_certificate(certificate)
end

function standardize_certificate(certificate)
    certificate = strip(string(certificate), ['\n', ' '])
    certificate
end

function compute_basis_validation_hash(result)
    ring, polys = parse_system_naive(result)
    certificate = get_certificate(ring, polys)
    standardize_certificate(certificate)
end

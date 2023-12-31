using LightXML
using AbstractAlgebra
using Groebner

function parse_polys_with_given_ring(ring, polys_str)
    polys_str = filter(!isempty, polys_str)
    base_field, vars = base_ring(ring), gens(ring)
    str_to_var_idx = Dict{String, Int}(string(v) => i for (i, v) in enumerate(vars))
    polys = Vector{AbstractAlgebra.elem_type(ring)}()
    constant_type = base_field == AbstractAlgebra.QQ ? Rational{BigInt} : Int64
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
            ring,
            constant_type,
            terms_exploded_str,
            str_to_var_idx
        )
        poly = ring(cfs, exps)
        push!(polys, poly)
    end
    polys
end

function load_system_symbolic_data(filename)
    # parse ex1.xml:
    # xdoc is an instance of XMLDocument, which maintains a tree structure
    xdoc = parse_file(filename)

    # get the root element
    xroot = LightXML.root(xdoc)  # an instance of XMLElement
    # print its name

    vars_str = []
    polys_str = []

    # traverse all its child nodes and print element names
    for c in child_nodes(xroot)  # c is an instance of XMLNode
        if is_elementnode(c)
            e = XMLElement(c)  # this makes an XMLElement instance
            if name(e) == "vars"
                vars_str = map(strip, split(LightXML.content(e), ","))
            elseif name(e) == "basis"
                for poly_node in child_nodes(e)
                    push!(polys_str, strip(LightXML.content(poly_node)))
                end
            end
        end
    end
    @assert !isempty(vars_str)
    @assert !isempty(polys_str)

    ring, xs = polynomial_ring(QQ, vars_str, ordering=:degrevlex)
    polys = parse_polys_with_given_ring(ring, polys_str)

    return ring, polys
end

ring, sys = load_system_symbolic_data(
    "/home/demin/SymbolicData/XMLResources/IntPS/SignalTheory.f966.xml"
)

@time groebner(sys);

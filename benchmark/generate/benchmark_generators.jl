
function generate_benchmark_source_for_groebner(name, system, dir, validate, nruns, time_filename)
    ring = parent(system[1])
    field = base_ring(ring)
    vars_repr = join(map(string, gens(ring)), ", ")
    field_repr = if iszero(characteristic(field))
        "QQ"
    else
        "GF($(characteristic(field)))"
    end
    vars_repr_quoted = map(s -> "$s", map(string, gens(ring)))
    ring_repr = """ring, ($vars_repr) = polynomial_ring(
        $field_repr, 
        $vars_repr_quoted, 
        internal_ordering=:degrevlex
    )"""
    system_repr = join(map(s -> "\t" * s, map(repr, system)), ",\n")
    system_repr = "system = [\n$system_repr\n]"
    buf = IOBuffer()
    println(buf, "# $name")
    println(buf, "#! format: off")
    println(buf, "using AbstractAlgebra, Groebner")
    println(buf, "")
    println(buf, ring_repr)
    println(buf, system_repr)
    String(take!(buf))
end

function generate_benchmark_source_for_singular(name, system, dir, validate, nruns, time_filename)
    generate_benchmark_source_for_groebner(name, system, dir, validate, nruns, time_filename)
end

function generate_benchmark_source_for_maplefgb(name, system, dir, validate, nruns, time_filename)
    ring = parent(system[1])
    field = base_ring(ring)
    buf = IOBuffer()
    println(buf, "# $name")
    println(buf, "with(Groebner):")
    println(buf, "with(PolynomialIdeals):")
    println(buf, "kernelopts(numcpus=1);")
    system_repr = replace(join(map(s -> "\t\t" * s, map(repr, system)), ",\n"), "//" => "/")
    vars_repr = join(map(string, gens(ring)), ", ")
    println(buf, "")
    println(buf, "runtime := 2^1000:")
    println(buf, "for i from 1 by 1 to 1 do")
    println(buf, "\tJ := [\n$system_repr\n\t]:")
    println(buf, "\tprint(\"Running $name\");")
    println(buf, "\tst := time[real]():")
    println(
        buf,
        "\tG := Groebner[Basis](J, tdeg($vars_repr), method=fgb, characteristic=$(characteristic(field))):"
    )
    println(buf, "\tprint(\"$name: \", time[real]() - st):")
    println(buf, "\truntime := min(runtime, time[real]() - st):")
    println(buf, "end do:")
    println(buf, "")
    println(buf, "timings_fn := \"$time_filename\":")
    println(buf, "FileTools[Text][WriteLine](timings_fn, \"$name\");")
    println(buf, "FileTools[Text][WriteLine](timings_fn, cat(\"total_time, \", String(runtime))):")
    if validate
        println(buf)
        output_fn = output_filename()
        println(buf, "output_fn := \"$dir/$output_fn\":")
        println(buf, "FileTools[Text][WriteLine](output_fn, \"$vars_repr\");")
        println(buf, "FileTools[Text][WriteLine](output_fn, \"$(characteristic(field))\");")
        println(
            buf,
            """
            for poly in G do
                FileTools[Text][WriteLine](output_fn, cat(String(poly), \",\")):
            end do:
            """
        )
    end
    String(take!(buf))
end

function generate_benchmark_source_for_mgb(name, system, dir, validate, nruns, time_filename)
    ring = parent(system[1])
    field = base_ring(ring)
    buf = IOBuffer()
    println(buf, "# $name")
    println(buf, "with(Groebner):")
    println(buf, "with(PolynomialIdeals):")
    println(buf, "kernelopts(numcpus=1);")
    system_repr = replace(join(map(s -> "\t\t" * s, map(repr, system)), ",\n"), "//" => "/")
    vars_repr = join(map(string, gens(ring)), ", ")
    println(buf, "")
    println(buf, "runtime := 2^1000:")
    println(buf, "for i from 1 by 1 to 1 do")
    println(buf, "\tJ := [\n$system_repr\n\t]:")
    println(buf, "\tprint(\"Running $name\");")
    println(buf, "\tst := time[real]():")
    println(buf, "\tG := libmgb:-gbasis($(characteristic(field)),0,[$vars_repr],J):")
    println(buf, "\tprint(\"$name: \", time[real]() - st):")
    println(buf, "\truntime := min(runtime, time[real]() - st):")
    println(buf, "end do:")
    println(buf, "")
    println(buf, "timings_fn := \"$time_filename\":")
    println(buf, "FileTools[Text][WriteLine](timings_fn, \"$name\");")
    println(buf, "FileTools[Text][WriteLine](timings_fn, cat(\"total_time, \", String(runtime))):")
    if validate
        println(buf)
        output_fn = output_filename()
        println(buf, "output_fn := \"$dir/$output_fn\":")
        println(buf, "FileTools[Text][WriteLine](output_fn, \"$vars_repr\");")
        println(buf, "FileTools[Text][WriteLine](output_fn, \"$(characteristic(field))\");")
        println(
            buf,
            """
            for poly in G do
                FileTools[Text][WriteLine](output_fn, cat(String(poly), \",\")):
            end do:
            """
        )
    end
    String(take!(buf))
end

function generate_benchmark_source_for_msolve(name, system, dir, validate, nruns, time_filename)
    ring = parent(system[1])
    field = base_ring(ring)
    vars_repr = join(map(string, gens(ring)), ", ")
    system_repr = replace(join(map(repr, system), ",\n"), "//" => "/")
    buf = IOBuffer()
    println(buf, "$vars_repr")
    println(buf, "$(characteristic(field))")
    println(buf, system_repr)
    String(take!(buf))
end

function generate_benchmark_source_for_openf4(name, system, dir, validate, nruns, time_filename)
    ring = parent(system[1])
    field = base_ring(ring)
    buf = IOBuffer()
    println(buf, "// $name")
    println(buf, "#include <iostream>")
    println(buf, "#include <fstream>")
    println(buf, "#include <string>")
    println(buf, "#include <vector>")
    println(buf, "#include <libopenf4.h>")
    println(
        buf,
        """
        using namespace std;

        int main (int argc, char **argv)
        {
        """
    )
    vars_repr = join(map(s -> "\tvariableName.push_back(\"$s\");", map(repr, gens(ring))), "\n")
    system_repr = join(map(s -> "\tpolynomialArray.emplace_back(\"$s\");", map(repr, system)), "\n")
    println(
        buf,
        """
        \tvector<string> polynomialArray;
        \tvector<string> variableName;

        $vars_repr
        $system_repr

        \tvector<string> basis = groebnerBasisF4($(characteristic(field)), $(length(gens(ring))), variableName, polynomialArray, 1, 0);

        \tstd::cout << \"The basis contains \" << basis.size() << \" elements.\" << std::endl;
        """
    )
    if validate
        vars_ = join(map(repr, gens(ring)), ", ")
        output_fn = output_filename()
        println(buf)
        println(
            buf,
            """
            \tofstream output;
            \toutput.open(\"$dir/$output_fn\");
            \toutput << \"$vars_\" << endl;
            \toutput << \"$(characteristic(field))\" << endl;
            \tfor (size_t i = 0; i < basis.size(); i++) {
            \t\toutput << basis[i] << \",\" << endl;
            \t}
            \toutput.close();
            """
        )
    end
    println(buf, "\treturn 0;")
    println(buf, "}")
    String(take!(buf))
end

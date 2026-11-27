"""Runtime-sized bridge to Poly Mojo's public polynomial API."""
from apn_mojo import Rational
from polynomials import Budget, QQ, ModInt, Modulus, NoContext, CoefficientField, SparsePoly, Term, Monomial, MonomialOrder, groebner, check_groebner_basis
from std.sys import argv
from std.time import perf_counter_ns


struct Reader(Movable):
    var text: String
    var position: Int

    def __init__(out self, var text: String):
        self.text = text^
        self.position = 0

    def line(mut self) raises -> String:
        var end = self.text.find("\n", self.position)
        if end < 0:
            raise Error("truncated benchmark input")
        var value = String(self.text[byte=self.position:end])
        self.position = end + 1
        return value^

    def count(mut self) raises -> Int:
        var value = Int(self.line())
        if value < 0:
            raise Error("negative benchmark count")
        return value


def run_case[C: CoefficientField](mut reader: Reader, ctx: C.Context, characteristic: String, directory: String, name: String, nruns: Int, certified: Bool, validate: Bool) raises:
    var work = Budget()
    var names = List[String]()
    var axes = List[Int64]()
    var nvars = reader.count()
    for i in range(nvars):
        names.append(reader.line())
        axes.append(Int64(i))
    var order = MonomialOrder(MonomialOrder.GREVLEX, axes^, work)
    var generators = List[SparsePoly[C]]()
    for _ in range(reader.count()):
        var terms = List[Term[C]]()
        for _ in range(reader.count()):
            var coefficient = C.zero(ctx)
            comptime if C == QQ:
                coefficient = rebind_var[C](QQ(Rational(reader.line())))
            else:
                coefficient = rebind_var[C](ModInt(UInt64(reader.count()), rebind[Modulus](ctx)))
            var pairs = List[Tuple[Int64, UInt32]]()
            for axis in range(nvars):
                var exponent = reader.count()
                if exponent > Int(UInt32.MAX):
                    raise Error("exponent exceeds UInt32")
                if exponent:
                    pairs.append((Int64(axis), UInt32(exponent)))
            terms.append(Term[C](Monomial(pairs^, work), coefficient^))
        generators.append(SparsePoly[C].from_terms(terms^, order, ctx, work))
    if reader.position != reader.text.byte_length():
        raise Error("trailing benchmark input")
    var best = Float64.MAX
    var best_check = Float64.MAX
    # One untimed warm-up, matching the upstream Groebner.jl backend.
    var samples = open(directory + "/samples.csv", "w")
    samples.write("run,compute_seconds,verify_seconds\n")
    for iteration in range(nruns + 1):
        var budget = Budget()
        var start = perf_counter_ns()
        var result = groebner(generators, order, budget, certify=certified, trace=certified)
        var elapsed = Float64(perf_counter_ns() - start) / 1e9
        var check_time = 0.0
        if certified:
            var checking = Budget()
            start = perf_counter_ns()
            if not result.trace:
                raise Error("certified basis did not retain a trace")
            if not result.trace.value().replay(generators, result.basis, checking):
                raise Error("basis trace verification failed")
            if not check_groebner_basis(generators, result, checking):
                raise Error("Groebner criterion verification failed")
            check_time = Float64(perf_counter_ns() - start) / 1e9
        if iteration:
            best = min(best, elapsed)
            best_check = min(best_check, check_time)
            samples.write(String(iteration, ",", elapsed, ",", check_time, "\n"))
        if iteration == nruns and validate:
            with open(directory + "/output", "w") as output:
                for i in range(len(names)):
                    if i:
                        output.write(", ")
                    output.write(names[i])
                output.write(String("\n", characteristic, "\n"))
                for i in range(len(result.basis)):
                    ref p = result.basis[i]
                    if i:
                        output.write(",\n")
                    if not len(p.terms):
                        output.write("0")
                    for j in range(len(p.terms)):
                        if j:
                            output.write(" + ")
                        ref term = p.terms[j]
                        comptime if C == QQ:
                            output.write(String(rebind[QQ](term.coefficient).value))
                        else:
                            output.write(String(rebind[ModInt](term.coefficient).value))
                        for pair in term.monomial.pairs:
                            output.write(String("*", names[Int(pair[0])], "^", pair[1]))
                output.write("\n")
    samples.close()
    with open(directory + "/timings", "w") as output:
        output.write(String(name, "\ntotal_time, ", best, "\n"))
        if certified:
            output.write(String("verification_time, ", best_check, "\n"))


def main() raises:
    var args = List[String]()
    for arg in argv():
        args.append(String(arg))
    if len(args) != 7:
        raise Error("usage: driver INPUT DIRECTORY NAME NRUNS plain|certified true|false")
    if args[5] != "plain" and args[5] != "certified":
        raise Error("invalid certificate mode")
    if args[6] != "true" and args[6] != "false":
        raise Error("invalid validation flag")
    var nruns = Int(args[4])
    if nruns < 1:
        raise Error("nruns must be positive")
    var reader = Reader(open(args[1], "r").read())
    if reader.line() != "POLYMOJO GROEBNER 1":
        raise Error("invalid interchange version")
    var characteristic = reader.line()
    var prime = UInt64(Int(characteristic))
    if prime == 0:
        run_case[QQ](reader, NoContext(), characteristic, args[2], args[3], nruns, args[5] == "certified", args[6] == "true")
    else:
        var work = Budget()
        run_case[ModInt](reader, Modulus(prime, work), characteristic, args[2], args[3], nruns, args[5] == "certified", args[6] == "true")

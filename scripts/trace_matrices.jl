using Revise, Groebner, AbstractAlgebra

# 1. Create the system
k = GF(2^30 + 3)
sys = sort(Groebner.Examples.katsuran(4, k=k), by=leading_monomial)

# 2. Learn the trace
trace, gb_truth = groebner_learn(sys)

# 3. Get the F4 matrices
matrices = Groebner.trace_export_matrices(trace.recorded_traces[UInt32])

# 4. Inspect the data structure
function make_matrix(sys, matrix)
    x = gens(parent(sys[1]))
    make_poly(idx, mult) = sys[idx] * prod(x .^ mult)
    reducers = map(make_poly, matrix.reducers.index, matrix.reducers.multiplier)
    to_be_reduced = map(make_poly, matrix.to_be_reduced.index, matrix.to_be_reduced.multiplier)
    (reducers=reducers, to_be_reduced=to_be_reduced)
end

matrix_1 = make_matrix(sys, matrices[1])
@info "" matrix_1.reducers matrix_1.to_be_reduced

# 5. Construct a basis by reducing the matrices
function follow_the_trace(sys, matrices)
    gb = sys
    for matrix in matrices
        matrix = make_matrix(gb, matrix)
        reduced = broadcast(AbstractAlgebra.normal_form, matrix.to_be_reduced, Ref(matrix.reducers))
        gb = vcat(gb, reduced)
    end
    gb
end

gb = follow_the_trace(sys, matrices)

@assert groebner(gb) == gb_truth

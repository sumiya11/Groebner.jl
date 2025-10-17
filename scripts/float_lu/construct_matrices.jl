using Revise, Groebner, Nemo, SparseArrays, LinearAlgebra

function make_matrix(sys, step)
    make_poly(idx, mult) = (sys[idx][1] .+ Ref(mult), sys[idx][2])
    reducers = map(make_poly, step.reducers.index, step.reducers.multiplier)
    to_be_reduced = map(make_poly, step.to_be_reduced.index, step.to_be_reduced.multiplier)
    (reducers=reducers, to_be_reduced=to_be_reduced)
end

function F4_to_coo(rows)
  T = eltype(rows[1][2])
  I, J, V = Vector{Int}(), Vector{Int}(), Vector{T}()
  for (i, row) in enumerate(rows)
    idx, val = rows[i]
    for j in 1:length(idx)
      push!(I, i)
      push!(J, idx[j])
      push!(V, val[j])
    end
  end
  I,J,V
end

function get_sparse_row(poly, monom_to_col, T)
  idx = [monom_to_col[m] for m in f4_monoms(poly)]
  val = map(T, f4_coeffs(poly))
  idx, val
end

function get_poly_from_sparse_row(row, col_to_monom)
  R = parent(first(values(col_to_monom)))
  k = base_ring(R)
  idx, val = row
  R(sum(col_to_monom[idx[i]] * k(val[i]) for i in 1:length(idx)))
end

function get_f4_from_sparse_row(row, col_to_monom)
  idx, val = row
  monoms = [col_to_monom[idx[i]] for i in 1:length(idx)]
  monoms, val
end

get_true_reducers(i) = Groebner.__DATA[(:red_up, i)]
get_true_to_be_reduced(i) = Groebner.__DATA[(:red_low, i)]
get_true_col_to_monom(i) = Groebner.__DATA[(:col_to_monom, i)]

function make_sparse_matrix(sys, step)
  i, step = step
  # @info "" step
  if isempty(step.reducers.index) && isempty(step.to_be_reduced.index)
    return []
  end
  reducers, to_be_reduced = make_matrix(sys, step)
  reducers = sort(reducers, by=f4_lm, lt=f4_drl, rev=true)
  # to_be_reduced = sort(to_be_reduced, by=f -> (leading_monomial(f), length(f)))
  @assert issorted(reducers, by=f4_lm, lt=f4_drl, rev=true)
  # @assert issorted(to_be_reduced, by=f -> (leading_monomial(f), length(f)))
  to_be_reduced = sort(to_be_reduced, by=f -> f4_lm(f), lt=f4_drl, rev=false)
  @assert allunique(f4_lm.(reducers))
  @debug "" i
  @debug "" reducers to_be_reduced
  @debug "" get_true_reducers(i)
  @debug "" get_true_to_be_reduced(i)
  @assert map(first, reducers) == map(first, get_true_reducers(i))
  @assert map(first, to_be_reduced) == map(first, get_true_to_be_reduced(i))
  cols = collect(union(reduce(union, map(f -> Set(f4_monoms(f)), reducers); init=Set()), reduce(union, map(f -> Set(f4_monoms(f)), to_be_reduced); init=Set())))
  lead_cols = Set(map(f -> f4_lm(f), reducers))
  function col_cmp(c1, c2)
    c1 in lead_cols && c2 in lead_cols && return f4_drl(c2,c1)
    !(c1 in lead_cols) && !(c2 in lead_cols) && return f4_drl(c2,c1)
    c1 in lead_cols
  end
  sort!(cols, lt=col_cmp)
  monom_to_col = Dict(cols .=> collect(1:length(cols)))
  col_to_monom = Dict(monom_to_col[k] => k for k in keys(monom_to_col))
  @debug "" col_to_monom
  @debug "" get_true_col_to_monom(i)
  @assert col_to_monom == get_true_col_to_monom(i)
  # @info "" reducers to_be_reduced cols lead_cols
  rows_up   = map(f -> get_sparse_row(f, monom_to_col, T_lu), reducers)
  rows_down = map(f -> get_sparse_row(f, monom_to_col, T_lu), to_be_reduced)
  # @info "" monom_to_col rows_up rows_down
  rows_all = vcat(rows_up,rows_down)
  I,J,V = F4_to_coo(rows_all)
  m, n, R = length(rows_all), length(cols), length(lead_cols)
  M = sparse(I,J,V,m,n)
  @debug "" m n R
  @debug "" M
  A1 = M[1:R,1:R]
  A2 = M[1:R, R+1:n]
  A3 = M[R+1:m, 1:R]
  A4 = M[R+1:m, R+1:n]
  @debug "" A1 A2 A3 A4
  @debug "" collect(A1) # cond(Array(A1), 2)
  @debug "" collect(diag(A1))
  if size(A1) == (0,0)
    C = A4
  else
    C = A4 - A3 * (A1 \ collect(A2))
  end
  @debug "" C
  if false
    C_lu = LinearAlgebra.lu(
      collect(C), RowNonZero(), allowsingular=true)
    @debug "" C_lu.p C_lu.P
    @debug "" C_lu issuccess(C_lu)
    C_ref = C_lu.U[invperm(C_lu.p), :]
    for i in 1:size(C_ref, 1)
      C_ref[i, :] /= C_ref[i, findfirst(!iszero, C_ref[i, :])]
    end
    @debug "" C_ref
    @assert C_lu.L * C_lu.U ≈ collect(C)[C_lu.p, :]
  else
    C_ref = Nemo.rref(Nemo.matrix(collect(C)))[2]
  end
  C_F4 = [(R.+findall(!iszero, C_ref[i, :]), nonzeros(sparse(C_ref[i, :]))) for i in 1:size(C_ref, 1)]
  # @info "" C_F4
  polys = map(row -> get_f4_from_sparse_row(row, col_to_monom), C_F4)
  sort!(polys, by=f4_lm, lt=f4_drl)
  # @info "" polys
  # polys = map(f -> f / leading_coefficient(f), polys)
  # @info "" polys
  return polys
end

f4_monoms(f) = f[1]
f4_coeffs(f) = f[2]
f4_lm(f)     = f[1][1]
f4_lc(f)     = f[2][1]
f4_len(f)    = length(f[1])

function f4_drl(e1::Vector{T}, e2::Vector{T}) where {T}
  m1 = Groebner.monom_construct_from_vector(Vector{Int}, e1)
  m2 = Groebner.monom_construct_from_vector(Vector{Int}, e2)
  Groebner.monom_isless(m1, m2, DegRevLex())
end

f4_to_aa(R, f) = sum(prod(gens(R) .^ mc[1]) * mc[2] for mc in zip(f4_monoms(f), f4_coeffs(f)))

function aa_to_f4(poly)
  collect(exponent_vectors(poly)), collect(coefficients(poly))
end

get_true_monoms(iteration) = map(first, Groebner.__DATA[iteration])

function follow_the_trace(sys, steps)
    @assert issorted(map(f4_lm, sys), lt=f4_drl)
    @assert all(f -> isone(f4_lc(f)), sys)
    
    gb = sys
    @assert map(f4_monoms, gb) == get_true_monoms(1)
  
    for (i, step) in enumerate(steps)
        # matrix = make_matrix(gb, step)
        # Can be computed via lin-alg
        # __reduced = broadcast(AbstractAlgebra.normal_form, matrix.to_be_reduced, Ref(matrix.reducers))
    
        reduced = make_sparse_matrix(gb, (i, step))
        @debug "" reduced
        gb = vcat(gb, reduced)
        
        # for f in map(f4_monoms, gb)
        #   println(f)
        # end
        # println("==true:")
        # for f in get_true_monoms(i+1)
        #   println(f)
        # end
        # @assert map(f4_monoms, gb) == get_true_monoms(i+1)
    end
    gb
end

to_zp(f, x) = f(numerator(x)) // f(denominator(x))

T_lu = Nemo.QQ

c    = Groebner.Examples.hexapod(k=Nemo.QQ)
# R, (x,y) = polynomial_ring(QQ, [:x, :y], internal_ordering=:degrevlex)
# c    = [x + y + 1, y*x^2 + y + 1, y^2*x^3 + y + 1]
c    = sort(c, by=leading_monomial)
c    = map(f -> f / leading_coefficient(f), c)
c_zp = map(f -> map_coefficients(c -> to_zp(GF(2^30+3),c), f), c)

t, gb1 = groebner_learn(c_zp);
groebner_apply!(t, c_zp);

steps = Groebner.trace_export_matrices(t.recorded_traces[UInt32])

c_f4 = aa_to_f4.(c)
gb2 = follow_the_trace(c_f4, steps)
gb2_rounded = map(f -> (f4_monoms(f), map(c -> rationalize(Float64(c), tol=1e-12), f4_coeffs(f))), gb2)

gb2_rounded_aa = f4_to_aa.(Ref(parent(c[1])), gb2_rounded)


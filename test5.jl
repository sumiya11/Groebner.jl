using AbstractAlgebra, Groebner, PrettyTables

groups, sys = Groebner.random_sys_bi_hom(
    GF(2^30 + 3), # ground
    (3,3),        # groups
    [(2,2),(2,2),(2,2),(2,2),(2,2),(2,2)] # bidegrees
)

# Set the stage
Groebner.SEMIGROUP_ON[] = false;
orig_system, new_system, varmap, relations, vars, tags, _, groups = Groebner.transform(sys, groups)
last_tag = tags[end]
hom_vars = vars[end-1:end]

# Affine
sys_affine = map(f -> evaluate(f, hom_vars,[parent(last_tag)(1),parent(last_tag)(1)]), orig_system)
@time gb_affine = Groebner.groebner(sys_affine; linalg=:deterministic, monoms=:dense)
summary_affine = deepcopy(Groebner.DATA); time_affine = deepcopy(Groebner.TIME)

# Semigroup
Groebner.SEMIGROUP_ON[] = false
@time gb = Groebner.groebner_semigroup(new_system, varmap, relations, groups);
summary = deepcopy(Groebner.DATA); time = deepcopy(Groebner.TIME)
@assert Groebner.isgroebner(gb)

# Smoke test
# new_sys_dehom = map(f -> evaluate(f, [last_tag],[parent(last_tag)(1)]), vcat(new_system, relations));
# gb_dehom = groebner(new_sys_dehom)
# gb_dehom_sub1 = map(f -> evaluate(f, [last_tag],[parent(last_tag)(1)]), gb_dehom);
# gb_sub1 = map(f -> evaluate(f, [last_tag],[parent(last_tag)(1)]), gb);
# @assert all(iszero, Groebner.normalform(gb_dehom_sub1, gb_sub1))
# @assert all(iszero, Groebner.normalform(gb_sub1, gb_dehom_sub1))

@info "" time_affine
@info "" time

for summary in (summary_affine, summary)
    push!(summary[:i], :TOTAL)
    push!(summary[:pairs], sum(summary[:pairs]))
    push!(summary[:matrix_size], reduce(.+, summary[:matrix_size]))
end

mykeys = [:i, :pairs, :degree, :matrix_size, :useful_rows]
pretty_table((; zip(mykeys, summary_affine[k] for k in mykeys)...); 
    header=mykeys, minimum_columns_width=8, title="affine",
    crop=:none)

pretty_table((; zip(mykeys, summary[k] for k in mykeys)...); 
    header=mykeys, minimum_columns_width=8, title="semigroup",
    crop=:none)

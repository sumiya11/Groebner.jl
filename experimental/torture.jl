include("../../../code/Julia/rur.jl")

include("../../../code/Data/Systems/reimer7.jl")

sys_z = convert_sys_to_sys_z(sys);
nn = 28
use_block = false

dm, Dq, sys_T, _vars, linform =
    prepare_system(sys_z, nn, AbstractAlgebra.parent(sys[1]), use_block);

pr = prevprime(2^nn - 1)
arithm = Groebner.ArithmeticZp(UInt64, UInt32, pr)
graph, t_learn, t_v, q, i_xw, t_xw, pr, gb_expvecs =
    learn_zdim_quo(sys_T, UInt32(pr), arithm, linform);
backup = deepcopy(graph)

for i in 1:400
    if i < 290
        continue
    end
    pr = Primes.prevprime(pr - 1)
    @info "" i pr
    expvecs, cfs_zz = extract_raw_data(sys_T)
    redflag, cfs_zp = reduce_mod_p(cfs_zz, UInt32(pr))
    success, gro = Groebner.groebner_applyX!(graph, cfs_zp, UInt32(pr))
    if !success
        print("\n*** bad prime for Gbasis detected ***\n")
        # The object may be corrupted after the failure. Revive it.
        graph = backup
        backup = deepcopy(backup)
        continue
    end
end

###

include("C:\\Users\\User\\Desktop\\RUR\\code\\Julia\\rur.jl")

QQ = AbstractAlgebra.QQ
polynomial_ring = AbstractAlgebra.polynomial_ring
include("C:\\Users\\User\\Desktop\\RUR\\code\\Data\\Systems\\reimer7.jl")

qq = zdim_parameterization(sys);

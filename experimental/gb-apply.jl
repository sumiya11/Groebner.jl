using BenchmarkTools
macro my_profview(ex)
    :((VSCodeServer.Profile).clear();
    VSCodeServer.Profile.init(n=10^8, delay=0.00001);
    VSCodeServer.Profile.start_timer();
    $ex;
    VSCodeServer.Profile.stop_timer();
    VSCodeServer.view_profile(;))
end

using AbstractAlgebra, Nemo

k = Groebner.henrion7(ground=AbstractAlgebra.GF(2^31 - 1), ordering=:degrevlex, np=Nemo)

context, gb = Groebner.groebner_learn(k);
context

@time flag, gb2 = Groebner.groebner_apply!(context, k);

gb == gb2

@my_profview for _ in 1:3
    Groebner.groebner_apply!(context, k)
end

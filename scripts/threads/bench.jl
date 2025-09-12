using Revise, Groebner, Nemo, TimerOutputs


@info "nthreads=$(Base.Threads.nthreads())"
sys = Groebner.Examples.alea6();

TimerOutputs.enable_timer!(Groebner._TIMER); reset_timer!(Groebner._TIMER); @time gb = groebner(sys; tasks=1); show(Groebner._TIMER, allocations=false); TimerOutputs.disable_timer!(Groebner._TIMER);

TimerOutputs.enable_timer!(Groebner._TIMER2); reset_timer!(Groebner._TIMER2); @time gb = groebner(sys; tasks=1); show(Groebner._TIMER2, allocations=false); println(); TimerOutputs.disable_timer!(Groebner._TIMER2);
TimerOutputs.enable_timer!(Groebner._TIMER2); reset_timer!(Groebner._TIMER2); @time gb = groebner(sys; tasks=2); show(Groebner._TIMER2, allocations=false); println(); TimerOutputs.disable_timer!(Groebner._TIMER2);
TimerOutputs.enable_timer!(Groebner._TIMER2); reset_timer!(Groebner._TIMER2); @time gb = groebner(sys; tasks=4); show(Groebner._TIMER2, allocations=false); println(); TimerOutputs.disable_timer!(Groebner._TIMER2);
TimerOutputs.enable_timer!(Groebner._TIMER2); reset_timer!(Groebner._TIMER2); @time gb = groebner(sys; tasks=8); show(Groebner._TIMER2, allocations=false); println(); TimerOutputs.disable_timer!(Groebner._TIMER2);

nothing

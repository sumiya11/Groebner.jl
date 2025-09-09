using Revise, Groebner, Nemo, TimerOutputs


@info "nthreads=$(Base.Threads.nthreads())"
sys = Groebner.Examples.chandran(10);

TimerOutputs.enable_timer!(Groebner._TIMER); reset_timer!(Groebner._TIMER); @time gb = groebner(sys; threaded=:no); show(Groebner._TIMER, allocations=false)

nothing
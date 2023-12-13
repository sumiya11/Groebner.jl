using NVTX, Groebner, AbstractAlgebra, Base.Threads

NVTX.@annotate function work()
    sleep(1)
end

NVTX.@range "A loop" begin
    Threads.@threads for i in 1:5
        work()
    end
end

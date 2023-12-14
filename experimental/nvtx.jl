using Base.Threads

function work()
    s = 0
    for i in 1:1000
        s += fetch(@spawn sum(rand(2^11)))
    end
    s
end

work()

@profview work()

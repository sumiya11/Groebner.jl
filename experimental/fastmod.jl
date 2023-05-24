
function normalsum(x, y, m)
    ((x % m) + (y % m)) % m
end

function fastsum(x, y, m)
    T = typeof(x)
    xy = x + y
    (xy - (Base.ashr_int(xy, sizeof(T) * 8 - 1) & (m^2))) % m
end

T = UInt8
m = T(11)
x = T(9)

for y in T(0):T(typemax(T) - x)
    normal = normalsum(x, y, m)
    fast = fastsum(x, y, m)
    if normal != fast
        println("$x + $y = $normal,  fast = $fast")
    end
end

################

T = UInt64
m = T(2^31 - 1)
x = T(232423)

for y in rand(T, 1000)
    normal = normalsum(x, y, m)
    fast = fastsum(x, y, m)
    if normal != fast
        println("$x + $y = $normal,  fast = $fast")
    end
end

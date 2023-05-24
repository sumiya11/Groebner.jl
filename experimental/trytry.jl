
function compute(x::UInt8)
    UInt8(2 * BigInt(x))
end

function compute(x::UInt64)
    UInt64(2 * BigInt(x))
end

function compute(x::BigInt)
    BigInt(2 * BigInt(x))
end

function owo(y::Integer)
    try
        x = UInt8(y)
        return Int(compute(x))
    catch E
        x = UInt64(y)
        return Int(compute(x))
    end
end

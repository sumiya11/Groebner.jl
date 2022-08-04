
struct ExponentVector1
    chunk1::UInt64
end

struct ExponentVector2
    chunk1::UInt64
    chunk2::UInt64
end

struct ExponentVector3
    chunk1::UInt64
    chunk2::UInt64
    chunk3::UInt64
end

function evsum(e1::ExponentVector1, 
                e2::ExponentVector1)
   ExponentVector1(e1.chunk1 + e2.chunk1)
end

function evsum(e1::ExponentVector2, 
                e2::ExponentVector2)
    ExponentVector2(e1.chunk1 + e2.chunk1, 
                    e1.chunk2 + e2.chunk2)
end

function evsum(e1::ExponentVector3, 
                e2::ExponentVector3)
    ExponentVector3(e1.chunk1 + e2.chunk1, 
                    e1.chunk2 + e2.chunk2,
                    e1.chunk3 + e2.chunk3)
end

e1 = ExponentVector1(1)
e2 = ExponentVector1(2)

@time evsum(e1,e2)

e1 = ExponentVector2(1, 1)
e2 = ExponentVector2(2, 3)

@time evsum(e1,e2)

e1 = ExponentVector3(1, 1, 1)
e2 = ExponentVector3(2, 3, 4)

@time evsum(e1,e2)

function revealify(
        m::T, 
        totaldegbits,
        variablebits) where {T}
    n = div(sizeof(m)*8 - totaldegbits, variablebits)
    bs = bitstring(m)
    tdb = totaldegbits
    s = [
        bs[tdb + 1 + (i-1)*variablebits:tdb + i*variablebits] 
        for i in 1:n
    ]
    [bs[1:totaldegbits], s...]
end

function extractify(
        m::T,
        totaldegbits,
        variablebits) where {T}
    n = div(sizeof(m)*8 - totaldegbits, variablebits)

end

function compactify(
        U::Type{T}, 
        monom::Vector{I},
        totaldegbits,
        variablebits) where {T, I}
    n = length(monom)
    @assert sizeof(T)*8 >= totaldegbits + variablebits*n
    indent = sizeof(T)*8 - totaldegbits
    x = zero(T)
    for (i, d) in enumerate(monom)
        x |= (d * 2^((i - 1)*variablebits))
        x += (d * 2^indent)
    end
    x
end

function monomprod(m1::T, m2::T) where {T}
    m1 + m2 
end

T = UInt64
m = UInt64[4, 0, 3]
totaldegbits = T(16)
variablebits = T(8)

m1 = UInt64[4, 0, 3]
m2 = UInt64[0, 1, 2]

i1 = compactify(T, m1, totaldegbits, variablebits)
ri1 = revealify(i1, totaldegbits, variablebits)
i2 = compactify(T, m2, totaldegbits, variablebits)
ri2 = revealify(i2, totaldegbits, variablebits)

ii = monomprod(i1, i2)
rii = revealify(ii, totaldegbits, variablebits)

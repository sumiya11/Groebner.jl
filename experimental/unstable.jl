
function foo(arr1, maxpairs::Int=typemax(Int), select_all::Bool=false) where {Deg}
    npairs::Int = length(arr1)
    if !select_all
        npairs = 10
    end
    npairs = min(npairs, maxpairs)
    @assert npairs > 0

    ps = arr1
    deg = arr1[1]

    sort!(arr1)

    if npairs > maxpairs
        navailable = npairs
        npairs = maxpairs
        lastlcm = arr1[npairs]
        while npairs < navailable && ps[npairs + 1] == lastlcm
            npairs += 1
        end
    end

    deg, npairs
end



#=
@testset "Chineese remainder theorem" begin

    modular_images(a, ms) = map(m -> mod(a, m), ms)

    ms = [2, 3, 5, 7]
    P = prod(ms)

    nums = [rand(1:P-1) for _ in 1:10]

    for a in nums
        rs = modular_images(a, ms)
        acrt = crt(rs, ms)
        @test acrt == a
    end

end
=#

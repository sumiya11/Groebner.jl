
using .GroebnerBases: groebner, noon3, noon4, noon5, katsura6, rootn,
                    is_groebner



using Logging
@testset "Groebner with fglm" begin
    noon = noon3(ground=QQ)

    gb = GroebnerBases.groebner(noon, zerodim=true, loglevel=Logging.Debug, linalg=:sparse)
    @test is_groebner(gb)

end

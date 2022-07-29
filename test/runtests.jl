using AstroEOM
using Test

@testset "AstroEOM.jl" begin

    u = [1.1, 0.0, 0.0, 0.0, 0.1, 0.0]
    p = [0.01214]
    t = 0.0
    du = zeros(6)
    eom_crtbp(du, u, p, t)

    @test isapprox(du, [0.0, 0.1, 0.0, -0.46406516020019184, 0.0, -0.0], rtol=1e-15)

end

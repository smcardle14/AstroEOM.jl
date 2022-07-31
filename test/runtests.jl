using AstroEOM
using Test
using BenchmarkTools
using LinearAlgebra

@testset "AstroEOM.jl" begin

    p = [0.01214]
    t = 0.0

    u = [1.1, 0.0, 0.0, 0.0, 0.1, 0.0]
    du = zeros(6)
    eom_crtbp!(du, u, p, t)
    @test isapprox(du, [0.0, 0.1, 0.0, -0.46406516020019184, 0.0, -0.0], rtol=1e-15)

    u = [1.1, 0.0, 0.0, 0.1]
    du = zeros(4)
    eom_crtbp_xy!(du, u, p, t)
    @test isapprox(du, [0.0, 0.1, -0.46406516020019184, 0.0], rtol=1e-15)

    u = vcat([1.1, 0.0, 0.0, 0.1], reshape(LinearAlgebra.I(4),16,1))
    du = zeros(20)
    eom_crtbp_stm_xy!(du, u, p, t)
    @test isapprox(du[1:4], [0.0, 0.1, -0.46406516020019184, 0.0], rtol=1e-15)

    # Benchmarking
    u = [1.1, 0.0, 0.0, 0.0, 0.1, 0.0]
    du = zeros(6)
    @btime eom_crtbp!($du, $u, $p, $t)

    u = [1.1, 0.0, 0.0, 0.1]
    du = zeros(4)
    @btime eom_crtbp_xy!($du, $u, $p, $t)

    u = vcat([1.1, 0.0, 0.0, 0.1], reshape(LinearAlgebra.I(4),16,1))
    du = zeros(20)
    @btime eom_crtbp_stm_xy!($du, $u, $p, $t)

end

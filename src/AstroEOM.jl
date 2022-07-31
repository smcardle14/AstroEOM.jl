module AstroEOM

using StaticArrays
using LinearAlgebra

function eom_crtbp!(du, u, p, t)

    x, y, z, vx, vy, vz = u
    μ = p[1]

    # Distance from primary
    x1 = x+μ
    r1 = sqrt(x1*x1+y*y+z*z)
    oneByR1_3 = 1.0/r1^3

    # Distance from secondary
    x2 = x+μ-1.0
    r2 = sqrt(x2*x2+y*y+z*z)
    oneByR2_3 = 1.0/r2^3

    du[1] = vx
    du[2] = vy
    du[3] = vz
    du[4] =  2.0*vy + x - (1.0-μ)*x1*oneByR1_3 - μ*x2*oneByR2_3
    du[5] = -2.0*vx + y - (1.0-μ)*y*oneByR1_3  - μ*y*oneByR2_3
    du[6] =              -(1.0-μ)*z*oneByR1_3  - μ*z*oneByR2_3

end

function eom_crtbp_xy!(du, u, p, t)

    x, y, vx, vy = u
    μ = p[1]

    # Distance from primary
    x1 = x+μ
    r1 = sqrt(x1*x1+y*y)
    oneByR1_3 = 1.0/r1^3

    # Distance from secondary
    x2 = x+μ-1.0
    r2 = sqrt(x2*x2+y*y)
    oneByR2_3 = 1.0/r2^3

    du[1] = vx
    du[2] = vy
    du[3] =  2.0*vy + x - (1.0-μ)*x1*oneByR1_3 - μ*x2*oneByR2_3
    du[4] = -2.0*vx + y - (1.0-μ)*y*oneByR1_3  - μ*y*oneByR2_3

end

function eom_crtbp_stm_xy!(du, u, p, t)

    x, y, vx, vy = u
    μ = p[1]

    # Distance from primary
    x1 = x+μ
    r1 = sqrt(x1*x1+y*y)

    # Distance from secondary
    x2 = x+μ-1
    r2 = sqrt(x2*x2+y*y)

    du[1] = vx
    du[2] = vy
    du[3] =  2*vy + x - (1-μ)*x1/r1^3 - μ*x2/r2^3
    du[4] = -2*vx + y - (1-μ)*y/r1^3  - μ*y/r2^3

    # Propagate STM
    G11 = 1 - μ/r2^3 + 3*μ*x2^2/r2^5 + (μ-1)/r1^3 - 3*x1^2*(μ-1)/r1^5
    G12 = 3*y*(r1^5*μ*x2 - r2^5*x1*(μ-1))/r1^5/r2^5
    G22 = 1 - μ/r2^3 + 3*y^2*μ/r2^5 + (μ-1)/r1^3 - 3*y^2*(μ-1)/r1^5

    A = SMatrix{4,4}(  0,   0,  G11, G12,
                       0,   0,  G12, G22,
                       1,   0,    0,  -2,
                       0,   1,    2,   0)

    STM = SMatrix{4,4}(@view u[5:20])
    dSTM = A * STM
    du[5:20] = dSTM

    return nothing

end

export eom_crtbp!
export eom_crtbp_xy!
export eom_crtbp_stm_xy!

end

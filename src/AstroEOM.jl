module AstroEOM

function eom_crtbp(du, u, p, t)

    x, y, z, vx, vy, vz = u
    μ = p[1]

    # Distance from primary
    x1 = x+μ
    r1 = sqrt(x1*x1+y*y+z*z)

    # Distance from secondary
    x2 = x+μ-1.0
    r2 = sqrt(x2*x2+y*y+z*z)

    du[1] = vx
    du[2] = vy
    du[3] = vz
    du[4] =  2.0*vy + x - (1.0-μ)*x1/r1^3 - μ*x2/r2^3
    du[5] = -2.0*vx + y - (1.0-μ)*y/r1^3  - μ*y/r2^3
    du[6] =              -(1.0-μ)*z/r1^3  - μ*z/r2^3

end

export eom_crtbp

end

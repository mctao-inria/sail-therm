function adjoint2idealsail(theta)

    tTheta   = tan(theta)
    β    = atan((-3 + sign(tTheta) * sqrt(9 + 8 * tTheta^2)) / tTheta / 4)

    return β
end

function control_ideal(x)
    r  = x[1:2]
    v = x[3:4]
    normr = sqrt( r[1]^2 + r[2]^2 + 1e-10)
    normv = sqrt( v[1]^2 + v[2]^2 + 1e-10)
    acos_arg = ( v[1] * r[1] + v[2] * r[2]) / normv / normr
    if acos_arg > 1
        acos_arg = 1
    end
    if acos_arg < -1
        acos_arg = -1
    end
    theta    = acos(acos_arg)
    tTheta   = tan(theta)
    
    β        = atan((-3 + sign(tTheta) * sqrt(9 + 8 * tTheta^2)) / tTheta / 4)
    return β
end

function srpsail2D(x, β, epsilon)
    # SRP of the ideal solar sail in 2D
    normr = sqrt( x[1]^2 + x[2]^2 + 1e-10)

    fsrp     = [ 2 * epsilon * cos(β)^3; 
                 2 * epsilon * sin(β) * cos(β)^2]
    fsrp     = fsrp / normr^2
    return fsrp
end
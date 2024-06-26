function adjoint2idealsail(theta)

    tTheta   = tan(theta)
    β    = atan((-3 + sign(tTheta) * sqrt(9 + 8 * tTheta^2)) / tTheta / 4)

    return β
end

function control_ideal(x)
    r  = x[1:3]
    v = x[4:6]
    normr = sqrt( r[1]^2 + r[2]^2 + r[3]^2 )
    normv = sqrt( v[1]^2 + v[2]^2 + v[3]^2 )
    acos_arg = ( v[1] * r[1] + v[2] * r[2] + v[3] * r[3] ) / normv / normr
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

function srpsail2D(x, β)
    # SRP of the ideal solar sail in 2D
    normr    = (x[1]^2 + x[2]^2)^(1/2)

    fsrp     = [ 2 * epsilon * cos(β)^3; 
                 2 * epsilon * sin(β) * cos(β)^2;
                 0]
    fsrp     = fsrp / normr^2
    return fsrp
end
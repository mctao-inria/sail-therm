function adjoint2idealsail(theta)

    tTheta   = tan(theta)
    β    = atan((-3 + sign(tTheta) * sqrt(9 + 8 * tTheta^2)) / tTheta / 4)

    return β
end

function control_ideal(x)
    r        = x[1:3]
    v        = x[4:6]

    acos_arg = dot(v / norm(v), r / norm(r))
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
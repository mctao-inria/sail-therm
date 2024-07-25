function theta2alpha(theta, b)
    # Local variables
    Nmax = 20

    # Decompose the vector b
    b1, b2, b3 = b[1], b[2], b[3]

    # Calculate the initial alpha
    ttheta = tan(theta)
    alpha = atan((-3.0 + sign(ttheta) * sqrt(9.0 + 8.0 * ttheta^2)) / (ttheta * 4.0))
    
    if alpha < 0.0
        alpha += pi
    end

    # Solve implicit function using Newton's method
    for ii in 1:Nmax
        cAlpha = cos(alpha)
        sAlpha = sin(alpha)

        aux1 = b1 + 3.0 * b2 * cAlpha^2 + 2.0 * b3 * cAlpha
        aux2 = b2 * cAlpha + b3
        aux3 = 2.0 * b2 * cAlpha + b3

        # Implicit function
        num = sAlpha * aux1
        den = cAlpha^2 * aux2 - sAlpha^2 * aux3

        tPsi = num / den
        psi = atan(tPsi)
        
        if psi < 0.0
            psi += pi
        end

        f = psi - theta

        # Derivative of the implicit function
        dnum = cAlpha * aux1 - sAlpha^2 * (6.0 * b2 * cAlpha + 2.0 * b3)
        dden = -cAlpha * sAlpha * (2.0 * aux2 + 2.0 * aux3 + b2 * cAlpha) + 2.0 * b2 * sAlpha^3
        df = (dnum * den - num * dden) / (den^2 + num^2)

        # Update alpha using Newton's method
        alpha -= f / df
    end

    return alpha
end


function control_non_ideal(x)
    r  = x[1:2]
    v = x[3:4]
    normr = rnorm(x[1:2])
    normv = rnorm(x[3:4])
    acos_arg = ( v[1] * r[1] + v[2] * r[2]) / normv / normr
    if acos_arg > 1
        acos_arg = 1
    end
    if acos_arg < -1
        acos_arg = -1
    end
    theta    = acos(acos_arg)
    #tTheta   = tan(theta)
    
    β        = theta2alpha(theta, b)
    u = [cos(β), sin(β)]
    return u
end

function srpsailnonideail2D(x, u, epsilon, b)
    # SRP of the ideal solar sail in 2D
    normr = rnorm(x[1:2])
    cb = u[1]
    sb = u[2]

    fsrp     = [(b[1] + b[2]) * cb^2 + b[3] * cb + b[1] * sb^2; 
                 b[2] * sb * cb + b[3] * sb]
    fsrp     = epsilon .* cb .*  fsrp ./ normr^2
    return fsrp
end
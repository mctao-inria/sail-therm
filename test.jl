using SciMLBase
using LinearAlgebra
using OrdinaryDiffEq
using Plots

mu = 1

function kepler!(dx, x, p, t)
    # Kepler equation
    #mu      = pars(1);
    #r = x[1:3]
    #v = x[4:6]
    
    #dv = - mu / norm(r)^3 * r
    #dx = [v; dv]
    normr = (x[1]^2 + x[2]^2 + x[3]^3)^(1/2)
    dx[1] = x[4]
    dx[2] = x[5]
    dx[3] = x[6]
    dx[4] = - mu / normr^3 * x[1]
    dx[5] = - mu / normr^3 * x[2]
    dx[6] = - mu / normr^3 * x[3]
end

x0 = [0.15000000000000002, 0.0, 0.0, 0.0, 3.4050261230332914, 3.4050261230344264e-6]
tf = 0.001
tspan = (0.0,tf)
prob = ODEProblem(kepler!,x0,tspan)
sol = solve(prob, DP5(), abstol = 1e-13, reltol = 1e-13, saveat = range(0., tf, length=1000))
plot(sol,vars=(1,2,3))
using SciMLBase
using LinearAlgebra
using OrdinaryDiffEq
using Plots

mu = 1

function kepler!(dx, x)
    # Kepler equation
    #mu      = pars(1);
    r = x[1:3]
    v = x[4:6]
    
    dv = - mu / norm(r)^3 * r
    dx = [v; dv]
end

x0 = [0.15000000000000002, 0.0, 0.0, 0.0, 3.4050261230332914, 3.4050261230344264e-6]
tf = 200
tspan = (0.0,tf)
prob = ODEProblem(kepler!,x0,tspan)
sol = solve(prob, DP5(), abstol = 1e-13, reltol = 1e-13, saveat = range(0., tf, length=1000))
plot(sol,vars=(1,2,3))
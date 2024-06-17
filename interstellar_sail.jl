import Pkg
Pkg.add("OptimalControl")
Pkg.add("LinearAlgebra")
Pkg.add("ForwardDiff")
Pkg.add("MINPACK")
Pkg.add("Plots")

using OptimalControl
using LinearAlgebra
using ForwardDiff
using MINPACK
using Plots

# Definition of the optical parameters 
rho     = 0.88 #0.88        # Specular reflection coefficient
s       = 1 #0.94        # Diffuse reflection coefficient 
eps_f   = 0.05        # Emissivity front coeff 0.05
eps_b   = 0.55        # Emissivity back coeff 0.55
Bf      = 0.79        # Front Lambertian coeff 0.79
Bb      = 0.55        # Back Lambertian coeff 0.55
eps     = 0 #(Bf * eps_f - Bb * eps_b) / (eps_b + eps_f) 
b       = [1 - rho * s, 2 * rho * s, Bf * rho * (1 - s) + (1 - rho) * eps];

# Solar and space constants
AU           = 1.495978707e11            # m Astronautical unit and Length unit HERE
LU           = AU                        # Length Unit = AU
mu           = 1
#mu           = 132712440042e9 / LU^3    % m^3/s^2
TU           = sqrt(AU^3/132712440042e9) # Time Unit
Cs           = 1367 * TU^3               # W m^-2 solar flux at 1AU
Tds          = 3                         # K temperature of the Deep Space
sigma        = 5.67051e-8 * TU^3         # Stefan-Boltzmann's constant [W / m^2 /K^4]
light_speed  = 299792458 / LU * TU       # Speed of light

# Sail parameters
Am           = 45 / LU^2                 # Area-to-mass ratio
Area         = 20 * 6 / LU^2             # Vane area * quantity of vanes
mass         = Area / Am                 # Mass of the satellite
temp0        = 293.15                    # Initial temperature of the satellite
Csrp         = Cs / light_speed
epsilon      = Am * Csrp

function kepl2cart(a, e, i, RAAN, omega, theta, mu)  # Copy of the matlab function
    # Initial conditions
    r_orb      = a * (1 - e^2) / (1 + e * cos(theta)) .* [cos(theta); sin(theta); 0]
    v_orb      = sqrt(mu / a / (1 - e^2)) .* [- sin(theta); e + cos(theta); 0]

    ARAAN      = [cos(RAAN) sin(RAAN) 0; - sin(RAAN) cos(RAAN) 0; 0 0 1]
    Ai         = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)]
    Aomega     = [cos(omega) sin(omega) 0; - sin(omega) cos(omega) 0; 0 0 1]

    Atot       = ARAAN' * Ai' * Aomega'

    r          = Atot * r_orb # [m]
    v          = Atot * v_orb # [m / s]
    return r, v
end

# Initial orbit data
rpsail        = 0.15                     # Periapsis distance elliptic orbit
a0            = (AU/LU + rpsail  * AU / LU ) /2
e0            = 1 - rpsail/a0
# OE0           = [0, 1e-6, 0, a0, e0, 0]
r0, v0        = kepl2cart(a0, e0, 1e-6, 0, 0, 0, mu)
x0            = [r0; v0]                 # Initial state
rE, vE        = kepl2cart(1, 0, 1e-6, 0, 0, 0, mu)
xEarth        = [rE; vE]
  
# Heat parameters for Kapton material
spec_heat    = 989 / LU^2 * TU^2                     # J/kg/K
heat_cap     = spec_heat * mass / LU^2 * TU^2        # J/K
Tlim         = 750                                   # K
temp_constr  = Tlim^4 - Tds^4

opt_constr   = (1 - rho)/(eps_f + eps_b)
heat_constr  = Cs/sigma


# Definition of the pars vector
pars         = [mu; epsilon; b']
pars0         = [mu; 0; b']

# Integration (MATLAB)
t0       = 0
tf       = 3600 * 24 * 30 * 12 * 3.5 / TU

function srpsail2D(x, β)
    # SRP of the ideal solar sail in 2D
    normr    = (x[1]^2 + x[2]^2)^(1/2)

    fsrp     = [ 2 * epsilon * cos(β)^3; 
                 2 * epsilon * sin(β) * cos(β)^2;
                 0]
    fsrp     = fsrp / normr^2
    return fsrp
end

function F0(x)
    # Kepler equation
    #mu      = pars(1);
    r       = x[1:3]
    v       = x[4:6]
    
    dvdt    = - mu / norm(r)^3 * r
    dxdt    = [v; dvdt]
    return dxdt
end

function F1(x, β)
    f           = atan(x[2], x[1])
    rot_matrix  = [cos(f) -sin(f) 0;
                    sin(f) cos(f) 0;
                    0 0 1 ]
    dvdt    = rot_matrix * srpsail2D(x, β)
    dxdt    = [0; 0; 0; dvdt]
    return dxdt
end

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

function readMatrixFromFile(filename)
    # Open the file for reading
    file = open(filename, "r")
    
    # Initialize an empty array to store the matrix rows
    matrix = []

    # Loop through each line in the file
    for line in eachline(file)
        if !isempty(line)
            row = parse.(Float64, split(line))
            push!(matrix, row)
        end
    end
    
    # Close the file
    close(file)
    
    # Convert the array to a matrix
    return hcat(matrix...)
end

@def ocp begin
    t ∈ [ t0, tf ], time
    x ∈ R⁶, state
    β ∈ R, control
    -π/2 ≤ β(t) ≤ π/2 
    x(t0) == x0
    ẋ(t) == F0(x(t)) + F1(x(t), β(t)) 
    # maximise energy -mu/r + 1/2 v^2
    -mu/sqrt(x[1](tf)^2 + x[2](tf)^2 + x[3](tf)^2) + 1/2 * (x[4](tf)^2 + x[5](tf)^2 + x[6](tf)^2) → max
    cos(β(t))/(x[1](t)^2 + x[2](t)^2 + x[3](t)^2) * opt_constr * heat_constr + temp_constr ≤ 0 # temperature constraint
end 

# initial guess as functions of time
filename = "matrix.txt"
x = readMatrixFromFile(filename)
println(size(x))

x = [1 2 3 4 5 6; 1 2 3 4 5 6; 1 2 3 4 5 6]

u = control_ideal(x)
initial_guess = (state=x, control=u, variable=0.05)

# Direct
nlp_sol = solve(ocp, init=initial_guess)

plot(nlp_sol, size=(600, 450))

nlp_sol
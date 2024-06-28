#import Pkg
#Pkg.add("OptimalControl")
#Pkg.add("LinearAlgebra")
#Pkg.add("ForwardDiff")
#Pkg.add("MINPACK")
#Pkg.add("Plots")
#Pkg.add("Interpolations")

using OptimalControl
#using LinearAlgebra
#using ForwardDiff
#using MINPACK
using Plots
using Interpolations

include("kepl2cart.jl")
include("control_ideal.jl")

# Definition of the optical parameters 
rho = 0.88 # Specular reflection coefficient
s = 1 #0.94 # Diffuse reflection coefficient 
eps_f = 0.05 # Emissivity front coeff 0.05
eps_b = 0.55 # Emissivity back coeff 0.55
Bf = 0.79 # Front Lambertian coeff 0.79
Bb = 0.55 # Back Lambertian coeff 0.55
eps = 0 #(Bf * eps_f - Bb * eps_b) / (eps_b + eps_f) 
b = [1 - rho * s, 2 * rho * s, Bf * rho * (1 - s) + (1 - rho) * eps];

# Solar and space constants
AU = 1.495978707e11            # m Astronautical unit and Length unit HERE
LU = AU                        # Length Unit = AU
mu = 1  #mu = 132712440042e9 / LU^3    % m^3/s^2
TU = sqrt(AU^3/132712440042e9) # Time Unit
Cs = 1367 * TU^3               # W m^-2 solar flux at 1AU
Tds = 3                         # K temperature of the Deep Space
sigma = 5.67051e-8 * TU^3         # Stefan-Boltzmann's constant [W / m^2 /K^4]
light_speed = 299792458 / LU * TU       # Speed of light

# Sail parameters
Am = 45 / LU^2                 # Area-to-mass ratio
Area = 20 * 6 / LU^2             # Vane area * quantity of vanes
mass = Area / Am                 # Mass of the satellite
temp0 = 293.15                    # Initial temperature of the satellite
Csrp = Cs / light_speed
epsilon = 1e-2 * Am * Csrp

# Initial orbit data
rpsail = 0.15                     # Periapsis distance elliptic orbit
a0 = (AU/LU + rpsail  * AU / LU ) /2
e0 = 1 - rpsail/a0
r0, v0 = kepl2cart(a0, e0, 1e-6, 0, 0, 0, mu)
x0 = [r0; v0]                 # Initial state
rE, vE = kepl2cart(1, 0, 1e-6, 0, 0, 0, mu)
xEarth = [rE; vE]
  
# Heat parameters for Kapton material
spec_heat = 989 / LU^2 * TU^2                     # J/kg/K
heat_cap = spec_heat * mass / LU^2 * TU^2        # J/K
Tlim = 700                                   # K
temp_constr = Tlim^4 - Tds^4

opt_constr = (1 - rho)/(eps_f + eps_b)
heat_constr = Cs/sigma

# Definition of the pars vector
pars = [mu; epsilon; b']
pars0 = [mu; 0; b']

# Read of the initial guess from Matlab
filename = "matrix.txt"
file = open(filename, "r")
file_content = read(file, String)
close(file)
lines = split(file_content, "\n")
matrix_data = []
for line in lines
    # Split each line by spaces and convert to float
    if !isempty(line)
        row = [parse(Float64, x) for x in split(line)]
        push!(matrix_data, row)
    end
end

# Integration (MATLAB)
t0 = 0
tf = 3600 * 24 * 30 * 12 * 3.5 / TU 

function F0(x)
    # Kepler equation
    normr = sqrt( x[1]^2 + x[2]^2 + x[3]^2 )
    r = x[1:3]
    v = x[4:6]
    dv = -mu / normr^3 * r
    dx = [v; dv]
    return dx
end

function F1(x, β)
    #f = atan(x[2], x[1])
    nx = sqrt( x[1]^2 + x[2]^2 + x[3]^2 )
    cf = x[1] / nx
    sf = x[2] / nx
    rot_matrix = [cf -sf  0
                  sf  cf  0
                  0   0   1]
    dvdt = rot_matrix * srpsail2D(x, β)
    dxdt = [0; 0; 0; dvdt]
    return dxdt
end





@def ocp begin
    t ∈ [ t0, tf ], time
    #x = (r₁, r₂, r₃, v₁, v₂, v₃) ∈ R⁶, state
    x ∈ R⁶, state
    β ∈ R, control
    -10 ≤ x₁(t) ≤ 10,   (1)
    -10 ≤ x₂(t) ≤ 10,    (2)
    -1 ≤ x₃(t) ≤ 1,      (3)
    -10 ≤ x₄(t) ≤ 10,    (4)
    -10 ≤ x₅(t) ≤ 10,    (5)
    -1 ≤ x₆(t) ≤ 1,      (6)
    -π/2 ≤ β(t) ≤ π/2 
    x(t0) == x0
    ẋ(t) == F0(x(t)) + F1(x(t), β(t)) 
    #cos(β(t)) / ( r₁(t)^2 + r₂(t)^2 + r₃(t)^2 ) * opt_constr * heat_constr + temp_constr ≤ 0
    cos(β(t)) / ( x₁(t)^2 + x₂(t)^2 + x₃(t)^2 ) * opt_constr * heat_constr + temp_constr ≤ 0
    -mu / sqrt( x₁(tf)^2 + x₂(tf)^2 + x₃(tf)^2 ) + 1/2 * ( x₄(tf)^2 + x₅(tf)^2 + x₆(tf)^2 ) → max
end

sol = solve(ocp)
plot_sol = Plots.plot(sol, size=(900, 1200))
savefig(plot_sol, "figures/plot_sol_without_initial_guess.pdf");

sol = sol.state.(sol.times)
Nsol = length(sol)
plot_traj2D = Plots.plot([ sol[i][1] for i ∈ 1:Nsol ], [ sol[i][2] for i ∈ 1:Nsol ], size=(600, 600), label="direct without initial guess")
savefig(plot_traj2D, "figures/plot_traj_without_initial_guess.pdf");

plot_traj_matlab = Plots.plot!(matrix_data[2], matrix_data[3], size=(600, 600), label="local-optimal")


###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################





# Interpolation of the initial guess
# plot_traj_matlab = Plots.plot(matrix_data[2], matrix_data[3], size=(600, 600))
t_inter = matrix_data[1]
N = length(t_inter)
x_inter = [ [ matrix_data[i][j] for i ∈ 2:7 ] for j ∈ 1:N ]
u_inter = control_ideal.(x_inter)

itp1 = LinearInterpolation(t_inter, [ x_inter[i][1] for i ∈ 1:N ])
itp2 = LinearInterpolation(t_inter, [ x_inter[i][2] for i ∈ 1:N ])
itp3 = LinearInterpolation(t_inter, [ x_inter[i][3] for i ∈ 1:N ])
itp4 = LinearInterpolation(t_inter, [ x_inter[i][4] for i ∈ 1:N ])
itp5 = LinearInterpolation(t_inter, [ x_inter[i][5] for i ∈ 1:N ])
itp6 = LinearInterpolation(t_inter, [ x_inter[i][6] for i ∈ 1:N ])
itp_u = LinearInterpolation(t_inter, u_inter)

## Figures
# Create individual plots
# plotx1 = Plots.plot(t_inter, itp1.(t_inter), grid = "off", framestyle = :box, legend = false)
# xlabel!("Time [0]")
# ylabel!("r₁ [-]")

# plotx2 = Plots.plot(t_inter, itp2.(t_inter), grid = "off", framestyle = :box, legend = false)
# xlabel!("Time [0]")
# ylabel!("r₂ [-]")

# plotx3 = Plots.plot(t_inter, itp3.(t_inter), grid = "off", framestyle = :box, legend = false)
# xlabel!("Time [0]")
# ylabel!("r₃ [-]")

# plotx4 = Plots.plot(t_inter, itp4.(t_inter), grid = "off", framestyle = :box, legend = false)
# xlabel!("Time [0]")
# ylabel!("v₁ [-]")

# plotx5 = Plots.plot(t_inter, itp5.(t_inter), grid = "off", framestyle = :box, legend = false)
# xlabel!("Time [0]")
# ylabel!("v₂ [-]")

# plotx6 = Plots.plot(t_inter, itp6.(t_inter), grid = "off", framestyle = :box, legend = false)
# xlabel!("Time [0]")
# ylabel!("v₃ [-]")

# Tile the plots in a 3x2 grid
# plotall = plot(plotx1, plotx2, plotx3, plotx4, plotx5, plotx6, layout = (3, 2), size = (1600, 1600));

# Display the tiled plot
# display(plotall);
# savefig(plotall, "figures/plotall.pdf");

# Control 
# Plots.plot(t_inter, itp_u.(t_inter), grid = "off", framestyle = :box, legend = false)
# xlabel!("Time [0]")
# ylabel!("u [-]")
# savefig(plotall, "figures/control.pdf");

x(t) = [itp1(t), itp2(t), itp3(t), itp4(t), itp5(t), itp6(t)]
β(t)  = itp_u(t)
# Initial guess
initial_guess = (state=x, control=β)

# Direct
sol = solve(ocp, init=initial_guess)

plot_sol = Plots.plot(sol, size=(900, 1200))
savefig(plot_sol, "figures/plot_sol_with_initial_guess.pdf");

x_sol = sol.state.(sol.times)
Nsize = length(x_sol)
plot_traj2D = Plots.plot([ x_sol[i][1] for i ∈ 1:Nsize ], [ x_sol[i][2] for i ∈ 1:Nsize ], size=(600, 600), label="direct with initial guess")
savefig(plot_traj2D, "figures/plot_traj_with_initial_guess.pdf");
plot_traj_matlab = Plots.plot!(matrix_data[2], matrix_data[3], size=(600, 600), label="local-optimal")


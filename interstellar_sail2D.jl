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
using JLD2

include("kepl2cart.jl")
include("control_ideal2D.jl")

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
epsilon = Am * Csrp

# Initial orbit data
# rpsail = 0.15                     # Periapsis distance elliptic orbit
rpsail = 0.15                     # Periapsis distance elliptic orbit
a0 = (AU/LU + rpsail  * AU / LU ) /2
e0 = 1 - rpsail/a0
r0, v0 = kepl2cart(a0, e0, 1e-6, 0, 0, 0, mu)
x0 = [r0[1:2]; v0[1:2]]                 # Initial state
rE, vE = kepl2cart(1, 0, 1e-6, 0, 0, 0, mu)
xEarth = [rE[1:2]; vE[1:2]]

# r0, v0 = kepl2cart(1, 0, 1e-6, 0, 0, 0, mu)
# x0 = [r0[1:2]; v0[1:2]]    

# Heat parameters for Kapton material
spec_heat = 989 / LU^2 * TU^2                     # J/kg/K
heat_cap = spec_heat * mass / LU^2 * TU^2        # J/K
Tlim = 800                                   # K
temp_constr = Tds^4 - Tlim^4 #Tlim^4 - Tds^4 

opt_constr = (1 - rho)/(eps_f + eps_b)
heat_constr = Cs/sigma

# Sun
rSun = 0.00465047 #AU

# Integration 
t0 = 0
tf = 3600 * 24 * 30 * 12 * 3.5 / TU 

# Definition of the pars vector
#pars = [mu; epsilon; b']
#pars0 = [mu; 0; b']

# Read of the initial guess from Matlab
# filename = "matrix.txt"
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

# Interpolation of the initial guess
# plot_traj_matlab = Plots.plot(matrix_data[2], matrix_data[3], size=(600, 600))
# t_inter = matrix_data[1]
# N = length(t_inter)
# x_inter = [ [ matrix_data[i][j] for i ∈ 2:7 ] for j ∈ 1:N ]
# u_inter = control_ideal.(x_inter)

# itp1 = LinearInterpolation(t_inter, [ x_inter[i][1] for i ∈ 1:N ])
# itp2 = LinearInterpolation(t_inter, [ x_inter[i][2] for i ∈ 1:N ])
# itp3 = LinearInterpolation(t_inter, [ x_inter[i][3] for i ∈ 1:N ])
# itp4 = LinearInterpolation(t_inter, [ x_inter[i][4] for i ∈ 1:N ])
# itp5 = LinearInterpolation(t_inter, [ x_inter[i][5] for i ∈ 1:N ])
# itp6 = LinearInterpolation(t_inter, [ x_inter[i][6] for i ∈ 1:N ])
# itp_u = LinearInterpolation(t_inter, u_inter)

# Integration from a random point x_init 
t_inter = matrix_data[1]
N = length(t_inter)
x_inter = [ [ matrix_data[i][j] for i ∈ 2:7 ] for j ∈ 1:N ]
u_inter = control_ideal.(x_inter)

N_init = 350
time_init = t_inter[N_init]
x0 = [x_inter[N_init][1:2]; x_inter[N_init][4:5]]

N_final = N

t0 = time_init
tf = t_inter[N_final] 

itp1 = LinearInterpolation(t_inter[N_init:N_final], [ x_inter[i][1] for i ∈ N_init:N_final ])
itp2 = LinearInterpolation(t_inter[N_init:N_final], [ x_inter[i][2] for i ∈ N_init:N_final ])
# itp3 = LinearInterpolation(t_inter[N_init:N_final], [ x_inter[i][3] for i ∈ N_init:N_final ])
itp4 = LinearInterpolation(t_inter[N_init:N_final], [ x_inter[i][4] for i ∈ N_init:N_final ])
itp5 = LinearInterpolation(t_inter[N_init:N_final], [ x_inter[i][5] for i ∈ N_init:N_final ])
# itp6 = LinearInterpolation(t_inter[N_init:N_final], [ x_inter[i][6] for i ∈ N_init:N_final ])
itp_u = LinearInterpolation(t_inter[N_init:N_final], u_inter[N_init:N_final])


function F0(x)
    # Kepler equation
    normr = sqrt( x[1]^2 + x[2]^2 + 1e-10)
    r = x[1:2]
    v = x[3:4]
    dv = -mu / normr^3 * r
    dx = [v; dv]
    return dx
end

function F1(x, β)
    #f = atan(x[2], x[1])
    normr = sqrt( x[1]^2 + x[2]^2 + 1e-10)
    cf = x[1] / normr
    sf = x[2] / normr
    rot_matrix = [cf -sf; sf  cf]
    dvdt = rot_matrix * srpsail2D(x, β, epsilon)
    dxdt = [0; 0; dvdt]
    return dxdt
end


function temperature(x, β)
    temp = (Tds^4 + cos(β) / ( x[1]^2 + x[2]^2 + 1e-10) * opt_constr * heat_constr)^(1/4)
    return temp
end

@def ocp begin
    t ∈ [ t0, tf ], time
    x ∈ R⁴, state
    β ∈ R, control
    -30 ≤ x₁(t) ≤ 30
    -30 ≤ x₂(t) ≤ 30
    -30 ≤ x₃(t) ≤ 30
    -30 ≤ x₄(t) ≤ 30
    -π/2 * 0.8 ≤ β(t) ≤ π/2 * 0.8
    x(t0) == x0
    ẋ(t) == F0(x(t)) + F1(x(t), β(t)) 
    sqrt(x₁(t)^2 + x₂(t)^2) ≥ 0.01
    #cos(β(t)) / ( x₁(t)^2 + x₂(t)^2 + 1e-3^2) * opt_constr * heat_constr + temp_constr ≤ 0
    -mu / sqrt( x₁(tf)^2 + x₂(tf)^2 + 1e-3^2) + 1/2 * ( x₃(tf)^2 + x₄(tf)^2 ) → max
end

function ocp_t0(N_0, N_f)
    global t0 = t_inter[N_0]
    global tf = t_inter[N_f]
    global x0 = [x_inter[N_0][1:2]; x_inter[N_0][4:5]]
    @def ocp begin
        t ∈ [ t0, tf ], time
        x ∈ R⁴, state
        β ∈ R, control
        -30 ≤ x₁(t) ≤ 30
        -30 ≤ x₂(t) ≤ 30
        -30 ≤ x₃(t) ≤ 30
        -30 ≤ x₄(t) ≤ 30
        -π/2 * 0.8 ≤ β(t) ≤ π/2 * 0.8
        x(t0) == x0
        ẋ(t) == F0(x(t)) + F1(x(t), β(t)) 
        sqrt(x₁(t)^2 + x₂(t)^2) ≥ 0.01
        # cos(β(t)) / ( x₁(t)^2 + x₂(t)^2 + 1e-3^2) * opt_constr * heat_constr + temp_constr ≤ 0
        -mu / sqrt( x₁(tf)^2 + x₂(tf)^2 + 1e-3^2) + 1/2 * ( x₃(tf)^2 + x₄(tf)^2 ) → max
    end
    return ocp
end

###########################################################################################################################################
#                           WITHOUT INITIAL GUESS
###########################################################################################################################################

# sol = solve(ocp, max_iter = 5000)#, grid_size = 100)
# sol = solve(ocp, init=sol, max_iter = 5000, grid_size = 100)#, grid_size = 100)
function fun_plot_sol(sol)
    x_sol = sol.state.(sol.times)
    Nsol = length(x_sol)
    plot_traj2D = Plots.plot([ x_sol[i][1] for i ∈ 1:Nsol ], [ x_sol[i][2] for i ∈ 1:Nsol ], size=(600, 600), label="direct without initial guess")
    Plots.plot!(matrix_data[2], matrix_data[3], size=(600, 600), label="local-optimal")
    scatter!([x_sol[1][1]], [x_sol[1][2]], label="beginning of the optimised arc" )
    scatter!([x_sol[end][1]], [x_sol[end][2]], label="end of the optimised arc" )
    plot!(rSun * cos.(range(0, 2*π, 50)), rSun * sin.(range(0, 2*π, 50)), label="Sun", color="yellow" )
    scatter!([0], [0], label="Sun", color="yellow" )
    return plot_traj2D
end

x(t) = [itp1(t), itp2(t), itp4(t), itp5(t)]
β(t)  = itp_u(t)

initial_guess = (state=x, control=β)
# initial_guess = sol_converged
# initial_guess = sol
# sol_converged = sol
# sol = solve(ocp, init=initial_guess)#, grid_size = 200)


# ocp = ocp_t0(250)
# initial_guess = sol

sol = solve(ocp, init=initial_guess, max_iter = 5000, grid_size = 500)
sol = solve(ocp, max_iter = 10000, grid_size = 200)

plot_sol = Plots.plot(sol, size=(900, 1200))
savefig(plot_sol, "figures/plot_sol_without_initial_guess.pdf");

x_sol = sol.state.(sol.times)
Nsol = length(x_sol)
# Nsol = 453
plot_traj2D = Plots.plot([ x_sol[i][1] for i ∈ 1:Nsol ], [ x_sol[i][2] for i ∈ 1:Nsol ], size=(600, 600), label="direct without initial guess")
plot_traj_matlab = Plots.plot!(matrix_data[2], matrix_data[3], size=(600, 600), label="local-optimal")
scatter!([x_sol[1][1]], [x_sol[1][2]], label="beginning of the optimised arc" )
scatter!([x_sol[end][1]], [x_sol[end][2]], label="end of the optimised arc" )
plot!(rSun * cos.(range(0, 2*π, 50)), rSun * sin.(range(0, 2*π, 50)), label="Sun", color="yellow" )
scatter!([0], [0], label="Sun", color="yellow" )
savefig(plot_traj2D, "figures/plot_traj_without_initial_guess.pdf");

sol_converged = sol
sol = sol_converged

# N_init = 250
# x0 = x_inter[N_init]
# t0 = t_inter[N_init]

# 270
init_loop = sol
t0_list = []
obj_list = []
sol_list = []
for Nt0_local = 215:-1:210 #Nt0_local=350:-10:200
    # Nt0_local = 1 
    Ntf_local = 499
    ocp_loop = ocp_t0(Nt0_local, Ntf_local)
    for Ngrid = 700:50:700
        global sol_loop = solve(ocp_loop, init=init_loop, max_iter = 3000, grid_size = Ngrid, print_level=0)
        if sol_loop.iterations == 5000
            println("Iterations exceeded while doing the continuation on the time")
            break
        end
        global init_loop = sol_loop
        println("Time: $(Nt0_local), Objective: $(sol_loop.objective), Iteration: $(sol_loop.iterations)")
        # println("Time: $(Ntf_local), Objective: $(sol_loop.objective), Iteration: $(sol_loop.iterations)")
        p = fun_plot_sol(sol_loop)
        display(p)
        push!(t0_list, t0)
        push!(obj_list, sol_loop.objective)
        push!(sol_list, sol_loop)
    end
end
sol = sol_list[end-1]
plot_sol = Plots.plot(sol, size=(900, 1200))
sol = sol_loop



for i = 1:length(sol_list)
    sol = sol_list[i]
    p = fun_plot_sol(sol)
    display(p)
    sleep(1)
end

save_object("sol2D_030.jld2", sol_list)
# sol_list = load("sol2D.jld2")


# sol = solve(ocp, init= sol_converged , max_iter = 5000, grid_size = 500)
# sol = solve(ocp, init= sol, max_iter = 5000, grid_size = 600)



β_sol = sol.control.(sol.times)
plot_temperature = Plots.plot(sol.times, temperature.(x_sol, β_sol), size=(600, 600), label="sail temperature")
plot!([0, sol.times[end]], [Tlim, Tlim], label="temperature limit")
savefig(plot_temperature, "figures/plot_temperature.pdf");

energy_sol = -mu./sqrt.([x_sol[i][1] for i ∈ 1:Nsol].^2 + [x_sol[i][2] for i ∈ 1:Nsol].^2 ) + 1/2 * ([x_sol[i][3] for i ∈ 1:Nsol].^2 + [x_sol[i][4] for i ∈ 1:Nsol].^2)
energy_local_optimal = -mu./sqrt.(matrix_data[2].^2 + matrix_data[3].^2 + matrix_data[4].^2) + 1/2 * (matrix_data[5].^2 + matrix_data[6].^2 + matrix_data[7].^2)

plot_energy = Plots.plot(sol.times, energy_sol, size=(600, 600), label="orbital energy")
plot!(matrix_data[1], energy_local_optimal, label="orbital energy, local-optimal")
savefig(plot_energy, "figures/plot_energy.pdf");

normr = sqrt.([ x_sol[i][1] for i ∈ 1:Nsol ].^2 + [ x_sol[i][2] for i ∈ 1:Nsol ].^2)
plot_normr = Plots.plot(sol.times, normr, size=(600, 600), label="distance fron the Sun")
plot!(matrix_data[1], sqrt.(matrix_data[2].^2 + matrix_data[3].^2), label="orbital energy, local-optimal")
savefig(plot_normr, "figures/plot_distance_from_sun.pdf");


plot_normr = Plots.plot(sol.times[215:225], normr[215:225], size=(600, 600), label="distance fron the Sun")


# xend = x_sol[end]
# -mu/(sqrt(xend[1]^2 + xend[2]^2)) + 1/2 * (xend[3]^2 + xend[4]^2)
a(1-e)

sol_converged = sol

using OptimalControl
using NLPModelsIpopt
using Plots
using Interpolations
using JLD2, JSON3

include("kepl2cart.jl")
include("control_ideal2D_cossin.jl")

function rnorm(x; ε=1e-9)
    return sqrt(sum(x.^2) + ε^2)
end

# Definition of the optical parameters 
rho = 0.88 # Specular reflection coefficient
eps_f = 0.05 # Emissivity front coeff 
eps_b = 0.55 # Emissivity back coeff 4

# Solar and space constants
LU = 1.495978707e11  # (m) Astronautical unit 
mu = 1  # m^3/s^2
TU = sqrt(LU^3/132712440042e9) # Time Unit
Cs = 1367 * TU^3  # W m^-2 solar flux at 1AU
Tds = 3  # K temperature of the Deep Space
sigma = 5.67051e-8 * TU^3  # Stefan-Boltzmann's constant [W / m^2 /K^4]
light_speed = 299792458 / LU * TU # Speed of light

# Sail parameters
Am = 45 / LU^2  # Area-to-mass ratio
Area = 20 * 6 / LU^2 # Vane area * quantity of vanes
mass = Area / Am  # Mass of the satellite
temp0 = 293.15  # Initial temperature of the satellite
Csrp = Cs / light_speed
epsilon = Am * Csrp

# Initial orbit data
rpsail = 0.15                     # Periapsis distance elliptic orbit
a0 = (1 + rpsail) /2
e0 = 1 - rpsail/a0
r0, v0 = kepl2cart(a0, e0, 1e-6, 0, 0, 0, mu)
x00 = [r0[1:2]; v0[1:2]]                 # Initial state

# Heat parameters for Kapton material
spec_heat = 989 / LU^2 * TU^2                     # J/kg/K
heat_cap = spec_heat * mass / LU^2 * TU^2        # J/K
Tlim = 600                                   # K
temp_constr = Tds^4 - Tlim^4 #Tlim^4 - Tds^4 

opt_constr = (1 - rho)/(eps_f + eps_b)
heat_constr = Cs/sigma

# Integration 
t0 = 0
tf = 3600 * 24 * 30 * 12 * 3.5 / TU 

# Read of the initial guess from Matlab
filename = "matrix.txt"
file = open(filename, "r")
file_content = read(file, String)
close(file)
lines = split(file_content, "\n")
matrix_data = []
for line in lines
    if !isempty(line)
        row = [parse(Float64, x) for x in split(line)]
        push!(matrix_data, row)
    end
end

# Integration from a random point x_init 
t_inter = matrix_data[1]#/ rnorm(x00[1:2])^2
N = length(t_inter)
x_inter = [ [ matrix_data[i][j] for i ∈ 2:7 ] for j ∈ 1:N ]
u_inter = control_ideal.(x_inter)

N_init = 490
time_init = t_inter[N_init]
x0 = [x_inter[N_init][1:2]; x_inter[N_init][4:5]]
# time_grid_nonuniform = t_inter[N_init:N]
N_final = N
# xf = [x_inter[N_final][1:2]; x_inter[N_final][4:5]]
# t_rescaled = t_inter ./ ([ x_inter[i][1] for i ∈ 1:N ].^2 + [ x_inter[i][2] for i ∈ 1:N ].^2)
# t0 = time_init #/ (x00[1]^2 + x00[2]^2)
# tf = t_inter[N_final] #/ (x00[1]^2 + x00[2]^2) #(xf[1]^2 + xf[2]^2)
# ([ x_inter[i][1] for i ∈ 1:N ].^2 + [ x_inter[i][2] for i ∈ 1:N ].^2)
# itp1 = linear_interpolation(t_inter[N_init:N_final], [ x_inter[i][1] for i ∈ N_init:N_final ], extrapolation_bc=Line())
# itp2 = linear_interpolation(t_inter[N_init:N_final], [ x_inter[i][2] for i ∈ N_init:N_final ], extrapolation_bc=Line())
# itp3 = linear_interpolation(t_inter[N_init:N_final], [ x_inter[i][4] for i ∈ N_init:N_final ], extrapolation_bc=Line())
# itp4 = linear_interpolation(t_inter[N_init:N_final], [ x_inter[i][5] for i ∈ N_init:N_final ], extrapolation_bc=Line())
# itp_u = linear_interpolation(t_inter[N_init:N_final], u_inter[N_init:N_final], extrapolation_bc=Line())

itp1 = linear_interpolation(t_inter, [ x_inter[i][1] for i ∈ 1:N_final ], extrapolation_bc=Line())
itp2 = linear_interpolation(t_inter, [ x_inter[i][2] for i ∈ 1:N_final ], extrapolation_bc=Line())
itp3 = linear_interpolation(t_inter, [ x_inter[i][4] for i ∈ 1:N_final ], extrapolation_bc=Line())
itp4 = linear_interpolation(t_inter, [ x_inter[i][5] for i ∈ 1:N_final ], extrapolation_bc=Line())
itp_u = linear_interpolation(t_inter, u_inter, extrapolation_bc=Line())


function F0(x)
    # Kepler equation
    normr = rnorm(x[1:2])
    r = x[1:2]
    v = x[3:4]
    dv = -mu / normr^3 * r
    dx = [v; dv] .* normr^2
    return dx
end

function F1(x, u)
    normr = rnorm(x[1:2])
    cf = x[1] / normr
    sf = x[2] / normr
    rot_matrix = [cf -sf; sf  cf]
    dvdt = rot_matrix * srpsail2D(x, u, epsilon)
    dxdt = [0; 0; dvdt] .* normr^2
    return dxdt
end

function temperature(x, u)
    normr = rnorm(x[1:2])
    temp = (Tds^4 + u[1] / ( normr ) * opt_constr * heat_constr)^(1/4)
    return temp
end

@def ocp begin
    t ∈ [ t0, tf ], time
    x ∈ R⁴, state
    u ∈ R², control
    -30 ≤ x₁(t) ≤ 30
    -30 ≤ x₂(t) ≤ 30
    -30 ≤ x₃(t) ≤ 30
    -30 ≤ x₄(t) ≤ 30
    cos(π/2 * 0.9) ≤ u₁(t) ≤ 1
    # cos(π/2 * 0.5) ≤ u₁(t) ≤ 1
    # sin(-π/2 * 0.8) ≤ u₂(t) ≤ sin(π/2 * 0.8)
    -1 ≤ u₂(t) ≤ 1
    u₁(t)^2 + u₂(t)^2 ≤ 1  #≤ 1
    x(t0) == x0
    ẋ(t) == F0(x(t)) + F1(x(t), u(t)) 
    (u₁(t)) / ( rnorm(x(t))) * opt_constr * heat_constr + temp_constr ≤ 0
    -mu / sqrt( rnorm(x(tf))) + 1/2 * ( x₃(tf)^2 + x₄(tf)^2 ) → max
end

function ocp_t0(N_0, N_f)
    global t0 = t_inter[N_0]
    global tf = t_inter[N_f]
    global x0 = [x_inter[N_0][1:2]; x_inter[N_0][4:5]]
    # global t0 = t_inter[N_0] #/ (x0[1]^2 + x0[2]^2)
    # global tf = t_inter[N_final] #/ (x0[1]^2 + x0[2]^2) 

    @def ocp begin
        t ∈ [ t0, tf ], time
        x = [r₁, r₂, v₁, v₂ ] ∈ R⁴, state
        u ∈ R², control
        -30 ≤ x₁(t) ≤ 30
        -30 ≤ x₂(t) ≤ 30
        -30 ≤ x₃(t) ≤ 30
        -30 ≤ x₄(t) ≤ 30
        cos(-π/2 * 0.9) ≤ u₁(t) ≤ 1
        # sin(-π/2 * 0.8) ≤ u₂(t) ≤ sin(π/2 * 0.8)
        -1 ≤ u₂(t) ≤ 1
        u₁(t)^2 + u₂(t)^2 ≤ 1
        # u₁(t)^2 + u₂(t)^2 == 1
        x(t0) == x0
        ẋ(t) == F0(x(t)) + F1(x(t), u(t)) 
        u₁(t) / ( rnorm(x(t))) * opt_constr * heat_constr + temp_constr ≤ 0
        -mu / sqrt( rnorm(x(tf))) + 1/2 * ( x₃(tf)^2 + x₄(tf)^2 ) → max
    end
    return ocp
end

function fun_plot_sol(sol)
    x_sol = sol.state.(sol.times)
    Nsol = length(x_sol)
    plot_traj2D = Plots.plot([ x_sol[i][1] for i ∈ 1:Nsol ], [ x_sol[i][2] for i ∈ 1:Nsol ], size=(600, 600), label="direct without initial guess")
    Plots.plot!(matrix_data[2], matrix_data[3], size=(600, 600), label="local-optimal")
    scatter!([x_sol[1][1]], [x_sol[1][2]], label="beginning of the optimised arc" )
    scatter!([x_sol[end][1]], [x_sol[end][2]], label="end of the optimised arc" )
    scatter!([0], [0], label="Sun", color="yellow" )
    return plot_traj2D
end

###########################################################################################################################################
#                           SIMPLE CODE
###########################################################################################################################################



x(t) = [itp1(t), itp2(t), itp3(t), itp4(t)]
u(t)  = itp_u(t)

initial_guess = (state=x, control=u)

# time_grid_non_uniform = []
# append!(time_grid_non_uniform, range(t0, t0 + (tf-t0)/2, 100)[1:end-1])
# append!(time_grid_non_uniform, range(t0 + (tf-t0)/2, tf, 200))

sol = solve(ocp; init=initial_guess, grid_size = 400, max_iter = 3000)
sol = solve(ocp; grid_size = 400, max_iter = 3000)

# sol = solve(ocp; init=initial_guess, time_grid = time_grid_nonuniform)
# sol = solve(ocp; init=initial_guess, time_grid = time_grid_non_uniform)
# time_grid_nonuniform

plot_sol = Plots.plot(sol, size=(900, 1200))
savefig(plot_sol, "figures/plot_sol.pdf");

x_sol = sol.state.(sol.times)
Nsol = length(x_sol)
# Nsol = 200
plot_traj2D = Plots.plot([ x_sol[i][1] for i ∈ 1:Nsol ], [ x_sol[i][2] for i ∈ 1:Nsol ], size=(600, 600), label="direct without initial guess", linewidth = 2, color = "blue", seriestype = :scatter)
savefig(plot_traj2D, "figures/plot_traj_dots.pdf");
plot_traj2D = Plots.plot([ x_sol[i][1] for i ∈ 1:Nsol ], [ x_sol[i][2] for i ∈ 1:Nsol ], size=(600, 600), label="direct without initial guess", linewidth = 2, color = "blue")#, seriestype = :scatter)
plot_traj_matlab = Plots.plot!(matrix_data[2], matrix_data[3], size=(600, 600), label="local-optimal", linewidth = 1, color = "red")
scatter!([x_sol[1][1]], [x_sol[1][2]], label="beginning of the optimised arc" )
scatter!([x_sol[end][1]], [x_sol[end][2]], label="end of the optimised arc" )
scatter!([0], [0], label="Sun", color="yellow" )
savefig(plot_traj2D, "figures/plot_traj.pdf");

u_sol = sol.control.(sol.times)

# t_rescaled = sol.times * ([ x_sol[i][1] for i ∈ 1:Nsol ].^2 + [ x_sol[i][2] for i ∈ 1:Nsol ].^2)
t_rescaled = sol.times * (x_sol[1][1].^2 + x_sol[1][2] ^2)
plot_beta = plot(t_rescaled, asind.([u_sol[i][2] for i ∈ 1:Nsol]), label="control angle", linewidth = 2, color = "blue")
ylabel!("[deg]")
savefig(plot_beta, "figures/plot_beta.pdf");
 

plot_temperature = Plots.plot(t_rescaled, temperature.(x_sol, u_sol), size=(600, 600), label="sail temperature", linewidth = 2, color = "blue")
plot!([t_rescaled[1], t_rescaled[end]], [Tlim, Tlim], label="temperature limit", linewidth = 2, color = "red")
ylabel!("[K]")
savefig(plot_temperature, "figures/plot_temperature.pdf");

energy_sol = -mu./sqrt.([x_sol[i][1] for i ∈ 1:Nsol].^2 + [x_sol[i][2] for i ∈ 1:Nsol].^2 ) + 1/2 * ([x_sol[i][3] for i ∈ 1:Nsol].^2 + [x_sol[i][4] for i ∈ 1:Nsol].^2)
energy_local_optimal = -mu./sqrt.(matrix_data[2].^2 + matrix_data[3].^2 + matrix_data[4].^2) + 1/2 * (matrix_data[5].^2 + matrix_data[6].^2 + matrix_data[7].^2)

plot_energy = Plots.plot(t_rescaled , energy_sol, size=(600, 600), label="orbital energy", linewidth = 2, color = "blue")
plot!(matrix_data[1], energy_local_optimal, label="orbital energy, local-optimal", linewidth = 1, color = "red")
savefig(plot_energy, "figures/plot_energy.pdf");

normr = sqrt.([ x_sol[i][1] for i ∈ 1:Nsol ].^2 + [ x_sol[i][2] for i ∈ 1:Nsol ].^2)
plot_normr = Plots.plot(sol.times, normr, size=(600, 600), label="distance from the Sun", linewidth = 2, color = "blue")
plot!(matrix_data[1], sqrt.(matrix_data[2].^2 + matrix_data[3].^2), label="distance from the Sun, local-optimal", linewidth = 1, color = "red")
# plot!([0, sol.times[end]], [0.01, 0.01], label="constraint")
savefig(plot_normr, "figures/plot_distance_from_sun.pdf");

# ylims!((0,0.4))
# xlims!((12,13))
# xlims!((21,22))

###########################################################################################################################################
#                           CONTINUATION ON T0
###########################################################################################################################################
# 142-340:0.056151306532663314
# 341(12)-639:0.016722408026755176 +-
# 640(17) - 838 :0.02349067839196195
# x > 0 : 12.65217391304348 - 12.735785953177256 (7) && 0 - 0.102873 (70)
# sol = load("run_17_07/solution_4")
# newtimegrid = sol.times
# newtimegrid = []
# append!(newtimegrid, time_grid_non_uniform[1:141])
# append!(newtimegrid, time_grid_non_uniform[142:5:342])
# append!(newtimegrid, time_grid_non_uniform[343:1:457])
# append!(newtimegrid, time_grid_non_uniform[458:20:838])
# init_loop = initial_guess
init_loop = sol
# init_loop = sol_last_converged
sol_list = []
# time_grid_nonuniform = sol.times
# time_grid_non_uniform = []
# append!(time_grid_non_uniform, range(t0, 12, 150)[1:end-1]) #14 200
# append!(time_grid_non_uniform, range(12, 19.5, 200)[1:end-1]) #17 400
# append!(time_grid_non_uniform, range(19.5, tf, 100))
# 18.701935 - 18.762195
# 309 - 313
# time_grid_non_uniform1 = []
# append!(time_grid_non_uniform1, time_grid_nonuniform[4:308])
# append!(time_grid_non_uniform1, range(18.701935,18.762195, 10))
# append!(time_grid_non_uniform1, time_grid_nonuniform[314:end])

for Nt0_local = 480:-10:450
    ocp_loop = ocp_t0(Nt0_local, N)
    # Ngrid = 500
    # time_grid_non_uniform = []
    # append!(time_grid_non_uniform, range(t0, 13, 150)[1:end-1]) #14 200
    # append!(time_grid_non_uniform, range(13, 18.5, 200)[1:end-1]) #17 400
    # append!(time_grid_non_uniform, range(18.5, tf, 100))
    # global time_grid_nonuniform = pushfirst!(time_grid_nonuniform, t0)
    # for Ngrid = 2000:10:2000 #1650
    Ngrid = 500
        global sol_loop = solve(ocp_loop, init=init_loop, grid_size = Ngrid, display = false, max_iter = 3000)
        # global sol_loop = solve(ocp_loop, time_grid = time_grid_nonuniform, init=init_loop, display = false, max_iter = 3000)
        # global sol_loop = solve(ocp_loop, time_grid = newtimegrid, init=init_loop, display = false, max_iter = 3000)
        # if sol_loop.iterations == 5000
    #     println("Iterations exceeded while doing the continuation on the time")
    #     break
    # end
        # x_sol = sol_loop.state.(sol.times)
        # Nsol = length(x_sol)
        # x1 = [ x_sol[i][1] for i ∈ 1:Nsol ]
        # iii = sol_loop.times[findall(i->(i>0), x1)]
        # println("The first time of positive x is: $(iii[1]), The last is: $(iii[end])")
        if sol_loop.objective < 2.5 && sol_loop.iterations < 3000
            global init_loop = sol_loop
            save(sol_loop, filename_prefix="solution_$(Nt0_local)")
            export_ocp_solution(sol_loop, filename_prefix="solution_$(Nt0_local)")
        end
        println("Time: $(Nt0_local), Objective: $(sol_loop.objective), Iteration: $(sol_loop.iterations)")
        # p = fun_plot_sol(sol_loop)
        # display(p)
    # end
    push!(sol_list, sol_loop)
end

sol = sol_list[end]
# sol_145 = sol
# sol_last_converged = sol

# x1 = [ x_sol[i][1] for i ∈ 1:Nsol ]
# findall(i->(i>0), x1)

# save_object("sol_12_07_ENTIRE.jld2", sol)
# sol_200 = sol
# sol_120 = sol
# sol_save = sol_300
# sol_150 = sol

sol_save = sol

# JLD save / load
# save(sol_save, filename_prefix="solution_145")
# sol = load("run_17_07/solution_4")
# println("Objective from loaded solution: ", sol_reloaded.objective)
# sol = load("sol_12_07_ENTIRE")

# JSON export / read
# export_ocp_solution(sol_save, filename_prefix="solution_145")
# sol_json = import_ocp_solution("my_solution")
# println("Objective from JSON discrete solution: ", sol_json.objective)

# fun_plot_sol(sol)
# x_sol = sol.state.(sol.times)
# Nsol = length(x_sol)
# energy_sol = -mu./sqrt.([x_sol[i][1] for i ∈ 1:Nsol].^2 + [x_sol[i][2] for i ∈ 1:Nsol].^2 ) + 1/2 * ([x_sol[i][3] for i ∈ 1:Nsol].^2 + [x_sol[i][4] for i ∈ 1:Nsol].^2)
# plot_energy = Plots.plot(sol.times, energy_sol, size=(600, 600), label="orbital energy")
# plot!(matrix_data[1][1:N], energy_local_optimal[1:N], label="orbital energy, local-optimal")

# sol_loaded = jld2.load("sol_12_07_ENTIRE.jld2")
# sol_loaded = sol_loaded["single_stored_object"]

# 1621
# 1:100 ! 0.178474:1.4921288944444446
# 101:900 1.5053981358024693: 12.107521980864197
# 901:1000 ! 12.120791222222225: 13.434446116666667
# 1001:1621     13.447715358024691:   21.674645
sol.times[1001:1621]
time_grid_refined = []
append!(time_grid_refined, range(0.178474, 1.4921288944444446, 3*100))
append!(time_grid_refined, range(1.5053981358024693, 12.107521980864197, 900-101+1))
append!(time_grid_refined, range(12.120791222222225, 13.434446116666667, 5*(1000-901+1)))
append!(time_grid_refined, range(13.447715358024691, 21.674645, (1621-1001+1)))
# Nbegin = 900
# Nsol = 1000
# plot_traj2D = Plots.plot([ x_sol[i][1] for i ∈ Nbegin:Nsol ], [ x_sol[i][2] for i ∈ Nbegin:Nsol ], size=(600, 600), label="direct without initial guess")
# plot_traj_matlab = Plots.plot!(matrix_data[2], matrix_data[3], size=(600, 600), label="local-optimal")
# scatter!([x_sol[1][1]], [x_sol[1][2]], label="beginning of the optimised arc" )
# scatter!([x_sol[end][1]], [x_sol[end][2]], label="end of the optimised arc" )
# scatter!([0], [0], label="Sun", color="yellow" )






plot_traj2D = Plots.plot([ x_sol[i][1] for i ∈ 1:Nsol ], [ x_sol[i][2] for i ∈ 1:Nsol ], size=(600, 600), label="direct without initial guess")
plot_traj_matlab = Plots.plot!(matrix_data[2], matrix_data[3], size=(600, 600), label="local-optimal")
scatter!([x_sol[1][1]], [x_sol[1][2]], label="beginning of the optimised arc" )
scatter!([x_sol[end][1]], [x_sol[end][2]], label="end of the optimised arc" )
scatter!([0], [0], label="Sun", color="yellow" )
ylims!((-1,1))
xlims!((-1,1))
savefig(plot_traj2D, "figures/plot_traj_zoom.pdf");





plot(t-Inter, itp1(t_inter))



x(t) = [itp1(t), itp2(t), itp3(t), itp4(t)]
u(t)  = itp_u(t)
t_inter
itp1(t_inter)
plot(itp1(t_inter), itp2(t_inter))


t_inter[1:130]
range(0, 0.534247, 130) # premier arc (début) 500
range(0.55369, 18.216503, 160) # 2 arc 200
range(18.249981, 21.674645, 210) # la fin 500
t1 = 0.534247
time_grid = []
push!(time_grid, range(0, 0.534247, 500), range(0.534247, 18.216503, 500)[2:end-1], range(18.216503, 21.674645, 500 ))







###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
#                           INDIRECT
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################


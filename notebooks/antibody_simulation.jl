using MKL 
using DifferentialEquations, LinearAlgebra, SparseArrays, Symbolics, Peaks;
using CSV,DataFrames,JSON
using StatsBase
using SparseArrays
include("utils.jl");

# Set parameters

## Parameters unchanged in simulations 
fixed_paras = JSON.parsefile("../data/fixed_parameters.json")
r_limit = fixed_paras["r_limit"]["value"]    # Maximum radial distance, unit um
N_r = fixed_paras["N_r"]["value"]       # Number of radial grid points
rgrid = range(0, stop=r_limit, length=N_r)
dr = step(rgrid)
mat_div = radial_laplacian_3d_Neumann(N_r, r_limit)
## smooth parameter for the init concentration field
epsilon = fixed_paras["epsilon"]["value"]   
## Radius of the nuclues in simulations
R0 = fixed_paras["R0"]["value"]   
 
## Parameters changed in simulations,
## Can choose difference figures for the parameters
# Load data from parameters.json
data = JSON.parsefile("../data/parameters.json")

figure_paras = data["Fig_2b1"]
D = figure_paras["Diffusion_Constant"]["D_a"]["value"]
k_on = figure_paras["Association_Rate"]["k_on"]["value"]
k_off = figure_paras["Dissociation_Rate"]["k_off"]["value"]
c_a_const = figure_paras["Antibody_Concentration"]["c_a"]["value"]
c_b_const = figure_paras["Epitope_Concentration"]["c_b"]["value"]

ps_temp = (D, k_on, k_off)

# Setup initial concentration vector
u_init = init_concentration_3D_spherical_smoothed(rgrid, epsilon, c_a_const, c_b_const, R0)

# Define ODEProblem
t_end = figure_paras["Simulation_Time"]["t_end"]["value"]
save_timestep = fixed_paras["save_timestep"]["value"]   

# Create ODE problem
ode_prob = ODEProblem(
    (du, u, p, t) -> diffuse_spherical_diff(du, u, p, t),
    u_init,
    (0.0, t_end),
    ps_temp
)

# Define the solver
solver = Rodas5P()

# Solve the ODE problem
diffuse_spherical_sol = solve(ode_prob, solver, saveat = save_timestep)
# Check results
println("Simulation completed.")
println(size(diffuse_spherical_sol))  # Should include time dimension

# Write results
df = DataFrame(time = diffuse_spherical_sol.t, u=diffuse_spherical_sol.u);
CSV.write("~/data/diffuse_spherical_sol.csv", df)




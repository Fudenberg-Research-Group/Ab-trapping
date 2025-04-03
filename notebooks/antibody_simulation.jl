
using MKL 
using DifferentialEquations, LinearAlgebra, SparseArrays, Symbolics, Peaks;
using CSV,DataFrames,JSON
using StatsBase
using SparseArrays
include("utils.jl");
import json 

# Set parameters

## Parameters unchanged in simulations 
r_limit = 75    # Maximum radial distance, unit um
N_r = 7500      # Number of radial grid points
rgrid = range(0, stop=r_limit, length=N_r)
dr = step(rgrid)
mat_div = radial_laplacian_3d_Neumann(N_r, r_limit)
## smooth parameter for the init concentration field
epsilon = 0.2
## Radius of the nuclues in simulations
R0 = 5


## Parameters changed in simulations,
## Can choose difference figures for the parameters
# Define diffusion parameters
# Load data from parameters.json
with open('parameters.json', 'r') as file:
    data = json.load(file)
    
figure_paras = data["Fig_2b"]
D = figure_paras["Diffusion_Constant"]["D_a"]["value"]
k_on = figure_paras["Association_Rate"]["k_on"]["value"]
k_off = figure_paras["Dissociation_Rate"]["k_off"]["value"]
c_a_const = figure_paras["Antibody_Concentration"]["c_a"]["value"]
c_b_const = figure_paras["Epitope_Concentration"]["c_b"]["value"]


ps_temp = (D, k_on, k_off)

# Setup initial concentration vector
u_init = init_concentration_3D_spherical_smoothed(rgrid, epsilon, c_a_const, c_b_const, R0)

# X = u_init
# dX = similar(X)

# Define ODEProblem
t_end = 3600
save_timestep = 10

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

df = DataFrame(time = diffuse_spherical_sol.t, u=diffuse_spherical_sol.u);
CSV.write("diffuse_spherical_sol.csv", df)




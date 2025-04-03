
using MKL 
using DifferentialEquations, LinearAlgebra, SparseArrays, Symbolics, Peaks;
using CSV,DataFrames,JSON
using StatsBase
using SparseArrays
include("utils.jl");


# Set parameters

## Parameters unchanged in simulations 
r_limit = 75    # Maximum radial distance, unit um
N_r = 7500      # Number of radial grid points
rgrid = range(0, stop=r_limit, length=N_r)
dr = step(rgrid)
mat_div = radial_laplacian_3d_Neumann(N_r, r_limit)
epsilon = 0.2
R0 = 5

## Parameters changed in simulations 
# Define diffusion parameters
D = 40
k_on = 1.6e9
k_off = 24
ps_temp = (D, k_on, k_off)
u1_const = 1.328e-8
u2_const = 6.34e-5


# Setup initial concentration vector
u_init = init_concentration_3D_spherical_smoothed(rgrid, epsilon, u1_const, u2_const, R0)

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




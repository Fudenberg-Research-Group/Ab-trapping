
using MKL 
using DifferentialEquations, LinearAlgebra, SparseArrays, Symbolics, Peaks;
using CSV,DataFrames,JSON
using StatsBase
using SparseArrays

function radial_laplacian_3d_Neumann(N_r, r_max)
    rd = range(0, stop=r_max, length=N_r)
    dr = step(rd)
    ir_indices = Int[]
    jc_indices = Int[]
    data_values = Float64[]

    for ir in 1:N_r
        r = rd[ir]
        center_idx = ir
        center_val = 0.0

        # Neumann boundary condition at r = 0 (ir == 1)
        if ir == 1
            push!(ir_indices, center_idx)
            push!(jc_indices, center_idx)
            push!(data_values, -1/(dr^2))   # Center contribution
 
            if N_r > 1
                push!(ir_indices, center_idx)
                push!(jc_indices, center_idx + 1)  # Right neighbor
                push!(data_values, 1/(dr^2))       # Contribution from T_1
            end
            
        # Internal points
        elseif ir > 1 && ir < N_r
            if ir > 1  # Left neighbor
                push!(ir_indices, center_idx)
                push!(jc_indices, ir - 1)
                push!(data_values, 1/(dr^2) - 1/(dr * r))  # Radial contribution
            end 

            if ir < N_r  # Right neighbor
                push!(ir_indices, center_idx)
                push!(jc_indices, ir + 1)
                push!(data_values, 1/(dr^2) + 1/(dr * r))  # Radial contribution
            end
            
            # Center position contribution for internal points
            center_val -= 2/(dr^2)  # Total contribution from neighbors
            push!(ir_indices, center_idx)
            push!(jc_indices, center_idx)
            push!(data_values, center_val)
        end
        
        # Neumann boundary condition at r = r_max (last point)
        if ir == N_r
            push!(ir_indices, center_idx)
            push!(jc_indices, center_idx)
            push!(data_values, -1/(dr^2))  # Center contribution

            if N_r > 1
                push!(ir_indices, center_idx)
                push!(jc_indices, center_idx - 1)  # Left neighbor
                push!(data_values, 1/(dr^2))  # Contribution from T_{N-1}
            end
        end
    end 
    
    return sparse(ir_indices, jc_indices, data_values, N_r, N_r)
end

function init_concentration_3D_spherical_smoothed(rgrid, epsilon, u1_const, u2_const, R0)
    # rgrid is a vector of radial distances
    N_r = length(rgrid)
    # Initialize the concentration field with zeros
    u = zeros(N_r, 3)  # Assuming u has three components (u1, u2, u3)
    for i in 1:N_r
        r = rgrid[i]
        u[i, 1] = u1_const/2 * (1 + tanh((r^2 - R0^2) / epsilon))  # Smooth transition for u1
        u[i, 2] = u2_const/2 * (1 - tanh((r^2 - R0^2) / epsilon))  # Smooth transition for u2
        u[i, 3] = 0  # u3 remains zero
    end
    u
end

# Define the ODE function for spherical diffusion
function diffuse_spherical_diff(du, u, p, t)
    D, kon, koff = p
    @inbounds du[:, 1] = D * mat_div * u[:, 1] + koff * u[:, 3] - kon * u[:, 1].*u[:, 2]
    @inbounds du[:, 2] = + koff * u[:, 3] - kon * u[:, 1].*u[:, 2]
    @inbounds du[:, 3] = - koff * u[:, 3] + kon * u[:, 1].*u[:, 2]
end

# Set parameters
r_limit = 75    # Maximum radial distance
N_r = 7500      # Number of radial grid points
rgrid = range(0, stop=r_limit, length=N_r)
dr = step(rgrid)
mat_div = radial_laplacian_3d_Neumann(N_r, r_limit)

# Define diffusion parameters
D = 40
k_on = 1.6e9
k_off = 24
ps_temp = (D, k_on, k_off)

epsilon = 0.2
R0 = 5
u1_const = 1.328e-8
u2_const = 6.34e-5
# Setup initial concentration vector
u_init = init_concentration_3D_spherical_smoothed(rgrid, epsilon, u1_const, u2_const, R0)

X = u_init
dX = similar(X)

# Define ODEProblem
t_end = 3600
savedt = 10

# Create ODE problem
ode_prob = ODEProblem(
    (du, u, p, t) -> diffuse_spherical_diff(du, u, p, t),
    X,
    (0.0, t_end),
    ps_temp
)

# Define the solver
solver = Rodas5P()

# Solve the ODE problem
diffuse_spherical_sol = solve(ode_prob, solver, saveat=savedt)
# Check results
println("Simulation completed.")
println(size(diffuse_spherical_sol))  # Should include time dimension

df = DataFrame(time = diffuse_spherical_sol.t, u=diffuse_spherical_sol.u);
CSV.write("diffuse_spherical_sol.csv", df)




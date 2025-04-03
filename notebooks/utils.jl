

# utils.jl

## Simulation 

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


## Visualization  


function extract_data(u_data::Vector{String}, time, indicator)
    # Step 1: Remove the square brackets and split the string by semicolons into rows
    rows = split(strip(u_data[time], ['[', ']']), ";")
    # Step 2: Extract the second column (3rd element in each row) and store as a vector of strings
    second_column = [split(row)[indicator] for row in rows]
    # Step 3: Convert the second column from string to Float64
    numbers = parse.(Float64, second_column)
    # Return the second column as numbers
    return numbers
end

function generate_radius_2D_new(side_length, shape)
    # Create a grid of size `shape x shape`
    y, x = [i for i in 1:shape, j in 1:shape], [j for i in 1:shape, j in 1:shape]
    # Define the center of the matrix
    center = (shape ÷ 2, shape ÷ 2)
    # Calculate the Euclidean distance from the center for each point
    radius_image = sqrt.((x .- center[2]).^2 + (y .- center[1]).^2)
    # Normalize the distance such that the maximum distance equals side_length / 2
    max_distance = sqrt(2) * (shape ÷ 2)  # max distance to the corner of the matrix
    scaled_radius_image = (radius_image / max_distance) * (side_length / 2) * sqrt(2)  
    return scaled_radius_image
end

function apply_function_to_radius(radius_image, func)
    # Apply the function to the radius values
    return func.(radius_image)
end



####

const tab10_colors = [
    RGB(0.121568, 0.466667, 0.705882),  # Blue
    RGB(1.000000, 0.498039, 0.054902),  # Orange
    RGB(0.172549, 0.627451, 0.172549),  # Green
    RGB(0.839216, 0.152941, 0.156863),  # Red
    RGB(0.580392, 0.403922, 0.741176),  # Purple
    RGB(0.549020, 0.337255, 0.294118),  # Brown
    RGB(0.890196, 0.466667, 0.760784),  # Pink
    RGB(0.498039, 0.498039, 0.498039),  # Gray
    RGB(0.737255, 0.741176, 0.133333),  # Olive
    RGB(0.090196, 0.745098, 0.811765)   # Cyan
]

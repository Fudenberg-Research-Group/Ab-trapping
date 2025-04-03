

# utils.jl

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

# This script will use the SmallBodyGravity package to plot the gravitational potential of the
# asteroid Bennu on the plane y=0 (in the body-fixed frame)
using SmallBodyGravity
using DelimitedFiles
using Plots
using ProgressMeter

# Constants and polyhedron processing
bennu_mass = 7.329e10
polyhedron_file = "geographos.obj"
G = 6.6743e-11

process_polyhedron(polyhedron_file, bennu_mass)

# Bennu's density for the used mass and shape 
# (!!! use the density from the generated `polyhedron_properties.txt` file !!!)
sigma = 1177.0535049082152

# Load all necessary data files into variables
centroid_edges = readdlm("centroid_edges.dat")
centroid_faces = readdlm("centroid_faces.dat")
e_e = readdlm("e_e.dat")
edges = readdlm("edges.dat")
faces = readdlm("faces.dat")
vertex = readdlm("vertex.dat")
n_f = readdlm("n_f.dat")
n_f_e = readdlm("n_f_e.dat")
n_fp_e = readdlm("n_fp_e.dat")
r_e_1 = readdlm("r_e_1.dat")
r_e_2 = readdlm("r_e_2.dat")
r_f_1 = readdlm("r_f_1.dat")
r_f_2 = readdlm("r_f_2.dat")
r_f_3 = readdlm("r_f_3.dat")

# Convert units from kilometers to meters
centroid_edges .*= 1000  # m
centroid_faces .*= 1000  # m
vertex .*= 1000  # m
r_f_1 .*= 1000  # m
r_f_2 .*= 1000  # m
r_f_3 .*= 1000  # m
r_e_1 .*= 1000  # m
r_e_2 .*= 1000  # m
e_e .*= 1000  # m

# Prepare parameter tuple
p = (centroid_edges, centroid_faces, e_e, edges, faces, vertex, n_f, n_f_e, n_fp_e, r_e_1, r_e_2, r_f_1, r_f_2, r_f_3, G, sigma)

# Set up grid for x and z
x_range = -600:10.0:600  # x-coordinates (in meters)
z_range = -600:10.0:600  # z-coordinates (in meters)

# Initialize array to store gravitational potential U
U_grid = zeros(length(z_range), length(x_range))

# Create a progress bar
total_iterations = length(x_range) * length(z_range)
progress = Progress(total_iterations, desc="Calculating gravitational potential...")

# Loop over x and z, keeping y = 0, and calculate U
for (i, z) in enumerate(z_range)
        for (j, x) in enumerate(x_range)
                local r_vec = [x, 0.0, z]  # Field point
                U_grid[i, j], Grav_Acceleration, _, _, P = polyhedron_model(p, r_vec)  # Compute gravitational potential
                if P == 1 # so that the potential inside the body is not plotted
                        U_grid[i,j] = NaN
                end
                # Update progress
                next!(progress)
        end
end

# Plot the results
p = heatmap(
x_range ./ 1e3, z_range ./ 1e3, U_grid;
xlabel = "x [km]", ylabel = "z [km]",
title = "Gravitational Potential U(x, z)",
colorbar_title = "U (m²/s²)",
aspect_ratio = :equal,
c = :viridis,
size = (1200, 1000)  # Set the size of the plot
)

# Save the plot
savefig(p, "gravitational_potential.png")  # Use savefig to save the plot created with specified size

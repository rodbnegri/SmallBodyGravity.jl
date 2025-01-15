# This script will use the SmallBodyGravity package to calculate the orbit of a spacecraft
# around the asteroid Bennu (in an inertial frame)


using DifferentialEquations # to integrate the trajectory 
using FileIO # to use the polyhedron in the PlotlyJS
using GeometryBasics
using SmallBodyGravity
using DelimitedFiles # needed for the readlm
using OrdinaryDiffEq
using Plots
using PlotlyJS # for making the plot
using LinearAlgebra # for using norm and other linear algebra functions

function EoM(dX, X, p, t)
	omega, p_polyhedron = p
	
	r_vec = X[1:3]
	v_vec = X[4:6]
	# Rotation matrix from inertial to body-fixed frame
	R_I_to_BF = [cos(omega * t) -sin(omega * t) 0; sin(omega * t) cos(omega * t) 0; 0 0 1]
	
	# position in the body-fixed frame
	r_vec_BF = R_I_to_BF * r_vec
	
	# calculate the acceleration
	_, Acceleration, _, _, P = polyhedron_model(p_polyhedron, r_vec_BF)
	
	if P != 0
		println("!!! The spacecraft collided !!!")
	end
	
	# Equations of Motion
	dX[1:3] = v_vec
	dX[4:6] = R_I_to_BF' * Acceleration
end

omega = 4.296057 / 3600.0 # rotation period [rad/s]
mass = 7.329e10
polyhedron_file = "bennu.obj"
G = 6.6743e-11

# integration time
tf = 2 * 24.0 * 3600.0

process_polyhedron(polyhedron_file, mass)

# Bennu density for the used mass and shape 
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
p_polyhedron = (centroid_edges, centroid_faces, e_e, edges, faces, vertex, n_f, n_f_e, n_fp_e, r_e_1, r_e_2, r_f_1, r_f_2, r_f_3, G, sigma)
p = (omega, p_polyhedron)

# Initial conditions 
r0 = [0.0, 100.0, 400.0]  # initial position in meters
v0 = [-0.08, 0.08, 0.0]     # initial velocity in m/s
X0 = vcat(r0, v0)

# Time span for integration
tspan = (0.0, tf)

println("Wait, the integration will take a long time (the spherical harmonics model that is being ported to Julia will solve this)")
# Solve the differential equations with Vern9 solver
prob = ODEProblem(EoM, X0, tspan, p)
sol = solve(prob, Vern9(), reltol=1e-8, abstol=1e-8)

# Extract trajectory
trajectory = sol.u |> collect

# Create a matrix with the first three elements of each vector
r_traj = hcat([v[1:3] for v in trajectory]...)  # Matrix where each column is the first three elements of a vector
# Calculate the magnitude of each vector (column) in r_traj
magnitudes = vec(norm.(eachcol(r_traj)))
# Find the largest magnitude
largest_magnitude = maximum(magnitudes) / 1e3

# Load the asteroid shape from the .obj file
mesh = load(polyhedron_file)

# Extract vertex positions and faces from the mesh
vertices_mesh = mesh.position
faces_mesh = mesh.faces

# Convert vertex positions to a 2D array for easier plotting
vertices_array = collect(Point3.(vertices_mesh))
vertices_x = [v[1] for v in vertices_array]
vertices_y = [v[2] for v in vertices_array]
vertices_z = [v[3] for v in vertices_array]

# Create the 3D plot with PlotlyJS only
plot_data = [
PlotlyJS.scatter3d(
x=r_traj[1, :] / 1e3, 
y=r_traj[2, :] / 1e3, 
z=r_traj[3, :] / 1e3, 
mode="lines", 
name="Trajectory",
line=attr(color="blue")
),
PlotlyJS.mesh3d(
x=vertices_x, 
y=vertices_y, 
z=vertices_z, 
color="gray", 
opacity=0.5, 
name="Asteroid"
)
]

# Set axis labels and layout
layout = Layout(
scene=attr(
xaxis=attr(
title="X (km)", 
showgrid=true, 
scaleanchor="y",  # Link the x-axis scale to the y-axis
range=[-largest_magnitude, largest_magnitude]  # Adjust to appropriate limits for your data
),
yaxis=attr(
title="Y (km)", 
showgrid=true, 
scaleanchor="x",  # Link the y-axis scale to the x-axis
range=[-largest_magnitude, largest_magnitude]  # Adjust to appropriate limits for your data
),
zaxis=attr(
title="Z (km)", 
showgrid=true,
scaleanchor="x",  # Link the z-axis scale to the x-axis
range=[-largest_magnitude, largest_magnitude]  # Adjust to appropriate limits for your data
),
aspectmode="cube"  # Ensures equal aspect ratio
),
title="Trajectory and Asteroid"
)

# Create the figure and display it
fig = PlotlyJS.plot(plot_data, layout)
display(fig)

# Save the plot as a PNG file
PlotlyJS.savefig(fig, "trajectory_plot.png")

println("Press Enter to continue...")
readline()


using SmallBodyGravity
using DelimitedFiles


bennuMass = 7.329e10
polyhedronFile = "bennu.obj"
G = 6.6743e-11

process_polyhedron("bennu.obj",bennuMass)

# Bennu's density for the used mass and shape 
# (!!! you should use the density found in the generated `polyhedron_properties.txt` file !!!)
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

# Small Body shape models are usually in kilometers (and in fact the process_polyhedron assumes so) to meters
        centroid_edges = centroid_edges * 1000 # m
        centroid_facets = centroid_faces * 1000 # m
        vertex = vertex * 1000 # m
        r_f_1 = r_f_1 * 1000 # m
        r_f_2 = r_f_2 * 1000 # m
        r_f_3 = r_f_3 * 1000 # m
        r_e_1 = r_e_1 * 1000 # m
        r_e_2 = r_e_2 * 1000 # m
        e_e = e_e * 1000 # m


p = centroid_edges, centroid_faces, e_e, edges, faces, vertex, n_f, n_f_e, n_fp_e, r_e_1, r_e_2, r_f_1, r_f_2, r_f_3, G, sigma

r_vec = [2e3; 0.0; 0.0]

polyhedron_model(p, r_vec)

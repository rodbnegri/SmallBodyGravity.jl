#
# This code was developed by Dr. Rodolfo Batista Negri at the National Institute for Space Research, Brazil.
#
# Polyhedron Gravity Model
# This function calculates the gravitational potential, force, and Laplacian for a polyhedron
# based on the edge and face-based summation methods, as outlined in:
# Scheeres, D.J., 2016. Orbital motion in strongly perturbed environments: applications to asteroid, comet and planetary satellite orbiters. Springer.
#
# The polyhedron is assumed to have triangular faces, and the function uses various geometric properties
# such as the centroid of edges and faces, the distance from the field point to edges and faces, and the normals
# to the edges and faces. The gravitational potential is computed using the closed-form solutions for constant density polyhedra.
#
# Function Parameters:
# - p: A tuple containing various polyhedron parameters, such as:
#     centroid_edges: Centroid of each edge
#     centroid_faces: Centroid of each face
#     e_e: Length of each edge
#     edges: Information about the edges of the polyhedron
#     face: Indicates the index of the vertices which form the face
#     vertex: The coordinates of each vertex
#     n_f: Normal vectors to the faces
#     n_f_e: Edge normals corresponding to each face
#     r_e_1, r_e_2: Vertices associated with each edge
#     r_f_1, r_f_2, r_f_3: Vertices associated with each face
#     G: Gravitational constant
#     sigma: Density of the polyhedron
# - r_vec: The vector representing the field point where the potential, force, and Laplacian are to be calculated.
#
# Returns:
# - U: Gravitational potential
# - F: Gravitational acceleration
# - Lapl: Laplacian of the potential

function polyhedron_model(p, r_vec)

    # Unpack parameters from input tuple `p`
    centroid_edges, centroid_faces, e_e, edges, face, vertex, n_f, n_f_e, n_fp_e, r_e_1, r_e_2, r_f_1, r_f_2, r_f_3, G, sigma = p

    # Calculate normal vectors for each edge used in the Edge summation
    n_f_Ee = n_f[Int.(edges[:, 3]), :]  # Normal vector to the face associated with each edge
    n_fp_Ee = n_f[Int.(edges[:, 4]), :]  # Normal vector to the other face associated with each edge

    # Calculate normal vectors for each face used in the face summation
    n_f_Ff = n_f

    # Compute the distance from each edge's and face's centroid to the field point
    re = centroid_edges - repeat(r_vec, size(centroid_edges, 1), 1)
    rf = centroid_faces - repeat(r_vec, size(centroid_faces, 1), 1)

    # Distance from each edge's vertex to the field point
    re1 = sqrt.(sum((r_e_1 - repeat(r_vec, size(centroid_edges, 1), 1)).^2, dims=2))
    re2 = sqrt.(sum((r_e_2 - repeat(r_vec, size(centroid_edges, 1), 1)).^2, dims=2))

    # Distance from each face's vertex to the field point
    rf1 = r_f_1 - repeat(r_vec, size(centroid_faces, 1), 1)
    rf2 = r_f_2 - repeat(r_vec, size(centroid_faces, 1), 1)
    rf3 = r_f_3 - repeat(r_vec, size(centroid_faces, 1), 1)
    
    # Norm for each face's vertices
    rf1_norm = sqrt.(sum(rf1.^2, dims=2))
    rf2_norm = sqrt.(sum(rf2.^2, dims=2))
    rf3_norm = sqrt.(sum(rf3.^2, dims=2))

    # Initialize variables for the summation over edges and faces
    sum_e_U = 0.0
    sum_e_F = 0.0
    sum_f_U = 0.0
    sum_f_F = 0.0
    Lapl = 0.0

    # Loop over edges and calculate the edge-based contributions to the potential and force
    for m = 1:size(centroid_edges, 1)
        Ee = n_f_Ee[m, :] * n_f_e[m, :]' + n_fp_Ee[m, :] * n_fp_e[m, :]

        # Logarithmic factor for edge length
        L_e = log.((re1[m] + re2[m] + e_e[m]) / (re1[m] + re2[m] - e_e[m]))

        sum_e_U += re[m, :]' * Ee * re[m, :] * L_e
        sum_e_F += Ee * re[m, :] * L_e

        if m <= size(centroid_faces, 1)
            # Sum over faces, inside edges' sum to save computational time
            Ff = n_f_Ff[m, :] * n_f_Ff[m, :]'

            # Calculate the solid angle for the face using the cross product
            omega_f = 2 * atan(dot(rf1[m, :], cross(rf2[m, :], rf3[m, :])), 
                               rf1_norm[m] * rf2_norm[m] * rf3_norm[m] + 
                               rf1_norm[m] * dot(rf2[m, :], rf3[m, :]) +
                               rf2_norm[m] * dot(rf3[m, :], rf1[m, :]) + 
                               rf3_norm[m] * dot(rf1[m, :], rf2[m, :]))

            sum_f_U += rf[m, :]' * Ff * rf[m, :] * omega_f
            sum_f_F += Ff * rf[m, :] * omega_f

            # Update the Laplacian term
            Lapl -= G * sigma * omega_f  # Laplacian of the potential

        end
    end

    # Compute the gravitational potential and acceleration
    U = G * sigma / 2 * (sum_e_U - sum_f_U)  # Potential
    F = -G * sigma * (sum_e_F - sum_f_F)  # Acceleration

    return U, F, Lapl  # Return the potential, acceleration, and Laplacian

end

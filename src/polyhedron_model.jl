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
#     sigma: Density of the small body
# - r_vec: The vector representing the field point where the potential, force, and Laplacian are to be calculated.
#
# Returns:
# - U: Gravitational potential
# - F: Gravitational acceleration
# - Lapl: Laplacian of the potential

function polyhedron_model(p, r_vec)
    
    centroid_edges, centroid_faces, e_e, edges, facets, vertex, n_f, n_f_e, n_fp_e, r_e_1, r_e_2, r_f_1, r_f_2, r_f_3, G, sigma = p
    
    # Calculate normal versors for each face used in the Edge summation
    n_f_Ee = n_f[Int.(edges[:,3]),:]
    n_fp_Ee = n_f[Int.(edges[:,4]),:]
    # Calculate normal versors for each face used in the Facets summation
    n_f_Ff = n_f;
    # Distance from each edge's/face's centroid to field point
    r_vec = reshape(r_vec, 1, :) # make sure r_vec is in proper dimension for using the repeat (this can be improved in later versions)
    re = centroid_edges - repeat(r_vec,size(centroid_edges,1),1);
    rf = centroid_faces - repeat(r_vec,size(centroid_faces,1),1);
    # Distance from each edge's vertex to field point
    re1 = sqrt.(sum((r_e_1 - repeat(r_vec,size(centroid_edges,1),1)).^2,dims=2));
    re2 = sqrt.(sum((r_e_2 - repeat(r_vec,size(centroid_edges,1),1)).^2,dims=2));
    # Distance from each face's vertex to field point
    rf1 = r_f_1 - repeat(r_vec,size(centroid_faces,1),1);
    rf2 = r_f_2 - repeat(r_vec,size(centroid_faces,1),1);
    rf3 = r_f_3 - repeat(r_vec,size(centroid_faces,1),1);
    rf1_norm = sqrt.( sum( rf1.^2 , dims=2 ) );
    rf2_norm = sqrt.( sum( rf2.^2 , dims=2 ) );
    rf3_norm = sqrt.( sum( rf3.^2 , dims=2 ) );
    
    sum_e_U = 0.;
    sum_e_F = 0.;
    sum_f_U = 0.;
    sum_f_F = 0.;
    Lapl = 0.;
    for m = 1:size(centroid_edges,1) # Sum over edges
        
        Ee = n_f_Ee[m,:]*n_f_e[m,:]' + n_fp_Ee[m,:]*n_fp_e[m,:]'
        
        L_e = log.( ( re1[m] + re2[m] + e_e[m] ) / ( re1[m] + re2[m] - e_e[m] ));
        
        sum_e_U = sum_e_U + re[m,:]'*Ee*re[m,:]*L_e;
        
        sum_e_F = sum_e_F .+ Ee*re[m,:]*L_e;
        
        
        if m <= size(centroid_faces,1) # Sum over facets (inside edges' sum to save computational time)
            
            Ff = n_f_Ff[m,:]*n_f_Ff[m,:]'
            
            omega_f =  2 * atan( dot( rf1[m,:] , cross( rf2[m,:] , rf3[m,:] ) ) ,
            rf1_norm[m] * rf2_norm[m] * rf3_norm[m] + rf1_norm[m] * dot( rf2[m,:] , rf3[m,:] ) +
            rf2_norm[m] * dot( rf3[m,:] , rf1[m,:] ) + rf3_norm[m] * dot( rf1[m,:] , rf2[m,:] ));
            
            sum_f_U = sum_f_U + rf[m,:]'*Ff*rf[m,:]*omega_f;
            
            sum_f_F = sum_f_F .+ Ff*rf[m,:]*omega_f;
            
            Lapl = Lapl - G * sigma * omega_f; # Laplacian
            
        end
        
    end
    U = G * sigma / 2 * ( sum_e_U - sum_f_U ); # Potential
    
    F = - G * sigma * ( sum_e_F - sum_f_F ); # Acceleration
    
    return U, F, Lapl
end

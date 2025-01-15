

#
# This code was developed by Dr. Rodolfo Batista Negri at the National Institute for Space Research, Brazil.
#
# Polyhedron Gravity Model
# This function calculates the gravitational potential, acceleration, Hessian and Laplacian for a polyhedron
# based on the edge and face-based summation methods, as outlined in:
# [1] Scheeres, D.J., 2016. Orbital motion in strongly perturbed environments: applications to asteroid, comet and planetary satellite orbiters. Springer.
# The original work used in the above reference is:
# [2]Werner, Robert A., and Daniel J. Scheeres. "Exterior gravitation of a polyhedron derived and compared with harmonic and mascon gravitation representations of asteroid 4769 Castalia." Celestial Mechanics and Dynamical Astronomy 65 (1996): 313-344.
# A good source for more info is:
# [3]Park, Ryan S., Robert A. Werner, and Shyam Bhaskaran. "Estimating small-body gravity field from shape model and navigation data." Journal of guidance, control, and dynamics 33.1 (2010): 212-221.
#
# The polyhedron is assumed to have triangular faces, and the function uses various geometric properties
# such as the centroid of edges and faces, the distance from the field point to edges and faces, and the normals
# to the edges and faces. The gravitational potential is computed using the closed-form solutions for constant density polyhedra.
#
# Function Parameters:
# - p: A tuple containing various polyhedron parameters, such as (their names are in accordance with reference [1]):
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
# - r_vec: The vector representing the field point where the potential, force, and Laplacian are to be calculated (note that this should be defined in the body-fixed frame).
#
# Returns:
# - U: Gravitational potential
# - A: Gravitational acceleration
# - H: Hessian matrix
# - L: Laplacian of the potential
# - P: if:
#          0: r_vec is outside the body
#          1: r_vec is inside the body
#          2: r_vec is on a face of the body
#          3: r_vec is on a edge or vertex
#
function polyhedron_model(p, r_vec)
  
  centroid_edges, centroid_faces, e_e, edges, facets, vertex, n_f, n_f_e, n_fp_e, r_e_1, r_e_2, r_f_1, r_f_2, r_f_3, G, sigma = p
  
  # Find normal unit vectors for each face used in the Edge summation
  n_f_Ee = n_f[Int.(edges[:,3]),:]
  n_fp_Ee = n_f[Int.(edges[:,4]),:]
  # Find normal unit vectors for each face used in the Faces summation
  n_f_Ff = n_f;
  # Distance from each edge's/face's centroid to field point
  r_vec = reshape(r_vec, 1, :) # make sure r_vec is in proper dimension for using the repeat (this can be improved in later versions)
  re = centroid_edges - repeat(r_vec,size(centroid_edges,1),1)
  rf = centroid_faces - repeat(r_vec,size(centroid_faces,1),1)
  # Distance from each edge's vertex to field point
  re1 = sqrt.(sum((r_e_1 - repeat(r_vec,size(centroid_edges,1),1)).^2,dims=2))
  re2 = sqrt.(sum((r_e_2 - repeat(r_vec,size(centroid_edges,1),1)).^2,dims=2))
  # Distance from each face's vertex to field point
  rf1 = r_f_1 - repeat(r_vec,size(centroid_faces,1),1)
  rf2 = r_f_2 - repeat(r_vec,size(centroid_faces,1),1)
  rf3 = r_f_3 - repeat(r_vec,size(centroid_faces,1),1)
  rf1_norm = sqrt.( sum( rf1.^2 , dims=2 ) )
  rf2_norm = sqrt.( sum( rf2.^2 , dims=2 ) )
  rf3_norm = sqrt.( sum( rf3.^2 , dims=2 ) )
  
  sum_e_H = zeros(3,3)
  sum_f_H = zeros(3,3)
  sum_e_U = 0.
  sum_e_A = zeros(3)
  sum_f_U = 0.
  sum_f_A = zeros(3) 
  L = 0.
  for m = 1:size(centroid_edges,1) # Sum over edges
    
    Ee = n_f_Ee[m,:]*n_f_e[m,:]' + n_fp_Ee[m,:]*n_fp_e[m,:]' # eq. 2.52 in the Reference
    
    L_e = log.( ( re1[m] + re2[m] + e_e[m] ) / ( re1[m] + re2[m] - e_e[m] )) # eq. 2.54 in the Reference 
    
    sum_e_U = sum_e_U + re[m,:]'*Ee*re[m,:]*L_e # summation for the edges in eq. 2.49
    
    sum_e_A = sum_e_A + Ee*re[m,:]*L_e # summation for the edges in eq. 2.49
    
    sum_e_H = sum_e_H + Ee * L_e 
    
    
    if m <= size(centroid_faces,1) # Sum over facets (inside edges' sum to save computational time)
      
      Ff = n_f_Ff[m,:]*n_f_Ff[m,:]'
      
      omega_f =  2 * atan( dot( rf1[m,:] , cross( rf2[m,:] , rf3[m,:] ) ) / (
      rf1_norm[m] * rf2_norm[m] * rf3_norm[m] + rf1_norm[m] * dot( rf2[m,:] , rf3[m,:] ) +
      rf2_norm[m] * dot( rf3[m,:] , rf1[m,:] ) + rf3_norm[m] * dot( rf1[m,:] , rf2[m,:] )))
      
      
      sum_f_U = sum_f_U + rf[m,:]'*Ff*rf[m,:]*omega_f
      
      sum_f_A = sum_f_A .+ Ff*rf[m,:]*omega_f
      
      sum_f_H = sum_f_H + Ff * omega_f
      
      L = L - G * sigma * omega_f # Laplacian
      
    end
    
  end
  U = G * sigma / 2 * ( sum_e_U - sum_f_U ) # Potential
  
  A = - G * sigma * ( sum_e_A - sum_f_A ) # Acceleration
  
  H = G * sigma * ( sum_e_H - sum_f_H )
  
  gamma = - L / G / sigma
  epsilon = 1e-8
  if isapprox(gamma, 4*pi, atol=epsilon)
    P = 1
  elseif isapprox(gamma, 0.0, atol=epsilon)
    P = 0
  elseif isapprox(gamma, 2*pi, atol=epsilon)
    P = 2
  else
    P = 3
  end
  
  return U, A, L, H, P
end

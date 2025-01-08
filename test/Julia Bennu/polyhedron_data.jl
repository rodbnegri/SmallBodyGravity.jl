# ------------------------------------------------------------------------
# Code Developed by Dr. Rodolfo Batista Negri
#
# Based on the work of  Werner and Scheeres [1]. This code should not be shared
# or distributed without prior consent from Dr. Negri (rodolfobnegri@yahoo.com.br).
#
# If you use this code in any manner, whether integrally, partially, or as
# inspiration, please cite references [2] and [3]. For instance:
# "... used the code developed by Negri [2,3], which applies the asteroid model 
# of Werner and Scheeres [1]."
#
# Please note that Dr. Negri does not assume any responsibility for the misuse 
# of this code. He is unable to provide assistance or support for its implementation.
#
# [1] Werner, R.A. and Scheeres, D.J., 1996. Exterior gravitation of a polyhedron 
# derived and compared with harmonic and mascon gravitation representations of 
# asteroid 4769 Castalia. Celestial Mechanics and Dynamical Astronomy, 65, pp.313-344.
#
# [2] Negri, R. B., “A Study of Dynamics, Guidance, Navigation, and Control 
# Applied to Asteroid Deflection," Ph.D. thesis, National Institute for Space Research,
# São José dos Campos, 2022.
#
# [3] Negri, R. B. and Prado, A.F. B. A., 2022. Autonomous and Robust Orbit-Keeping 
# for Small-Body Missions. Journal of Guidance, Control, and Dynamics, 45(3), pp.587-598.
# ------------------------------------------------------------------------

function polyhedron_data(file_in::String, file_edge::String)
	
	
	# coords da edge, duas primeiras colunas dão os vertex de cada edge,
	# duas colunas restantes dão as faces que formam a edge
	edges = readdlm(file_edge,Int)
	
	faces = readdlm("faces.dat", Int)
	vertex = readdlm("vertex.dat")
	
	# Coordinates of each vertex of each face
	r_f_1 = vertex[faces[:,1], :]
	r_f_2 = vertex[faces[:,2], :]
	r_f_3 = vertex[faces[:,3], :]
	
	# Coordinates of each vertex of each edge
	r_e_1 = vertex[edges[:,1], :]
	r_e_2 = vertex[edges[:,2], :]
	
	# Size of each edge
	e_e_vec = r_e_2 .- r_e_1
	e_e = sqrt.(sum((r_e_2 .- r_e_1).^2, dims=2))
	
	# Corresponding faces for each edge
	f_e = faces[edges[:,3], :]
	fp_e = faces[edges[:,4], :]
	
	# Find centroid of each face and edge (it could be any point on both)
	centroid_faces = (r_f_1 .+ r_f_2 .+ r_f_3) ./ 3
	centroid_edges = (r_e_1 .+ r_e_2) ./ 2
	
	
	
	# Normal vector to each face
	n_f = zeros(size(r_f_1))
	# normal to each edge
	n_f_e = n_f
	n_fp_e = n_f
  #TODO: tem ALGUM ERROR AQUI
	for ii = 1:size(r_f_1,1)
		nonunitary_n_f = cross(r_f_2[ii,:] - r_f_1[ii,:], r_f_3[ii,:]- r_f_2[ii,:])
		n_f[ii,:] = nonunitary_n_f / norm(nonunitary_n_f)
		norm_btw_centroids = norm( centroid_edges[ii,:] - centroid_faces[edges[ii,3],:] )
		nn_paralel = (centroid_edges[ii,:] - centroid_faces[edges[ii,3],:]) / norm_btw_centroids# versor from centroid of the face to the edge's centroid
		nn_perpendicular = cross( nn_paralel , e_e_vec[ii,:])# vector perpendicular to edge and face
		n_f_e[ii,:] = cross(nn_perpendicular,e_e_vec[ii,:]) / norm(cross(nn_perpendicular,e_e_vec[ii,:]))  # versor parallel to face and perpendicular to edge
		theta = acosd(dot(nn_paralel , n_f_e[ii,:]))
    @show n_f_e[ii,:], theta
		if theta>90.0
			n_f_e = -n_f_e # making sure n_f_e is pointing outside its face
		end
    @show n_f_e[ii,:]
    sleep(10)

		# doing the same for n_fp_e
		norm_btw_centroids = norm( centroid_edges[ii,:] - centroid_faces[edges[ii,4],:] )
		nn_paralel = (centroid_edges[ii,:] - centroid_faces[edges[ii,4],:]) / norm_btw_centroids# versor from centroid of the face to the edge's centroid
		nn_perpendicular = cross( nn_paralel , e_e_vec[ii,:])# vector perpendicular to edge and face
		n_fp_e[ii,:] = cross(nn_perpendicular,e_e_vec[ii,:]) / norm(cross(nn_perpendicular,e_e_vec[ii,:]))  # versor parallel to face and perpendicular to edge
		theta = acosd(dot(nn_paralel , n_fp_e[ii,:]))
		if theta>90.0
			n_fp_e= -n_fp_e # making sure n_fp_e is pointing outside its face
		end
		
	end
	return	
	
	# Create data for other programs
	open("edges.dat", "w") do f_edges
		for row in eachrow(edges)
			@printf(f_edges, "%d\t%d\t%d\t%d\n", row...)
		end
	end
	
	open("r_f_1.dat", "w") do f_r_f_1
		for row in eachrow(r_f_1)
			@printf(f_r_f_1, "%.12e\t%.12e\t%.12e\n", row...)
		end
	end
	
	open("r_f_2.dat", "w") do f_r_f_2
		for row in eachrow(r_f_2)
			@printf(f_r_f_2, "%.12e\t%.12e\t%.12e\n", row...)
		end
	end
	
	open("r_f_3.dat", "w") do f_r_f_3
		for row in eachrow(r_f_3)
			@printf(f_r_f_3, "%.12e\t%.12e\t%.12e\n", row...)
		end
	end
	
	open("r_e_1.dat", "w") do f_r_e_1
		for row in eachrow(r_e_1)
			@printf(f_r_e_1, "%.12e\t%.12e\t%.12e\n", row...)
		end
	end
	
	open("r_e_2.dat", "w") do f_r_e_2
		for row in eachrow(r_e_2)
			@printf(f_r_e_2, "%.12e\t%.12e\t%.12e\n", row...)
		end
	end
	
	open("e_e.dat", "w") do f_e_e
		for value in e_e
			@printf(f_e_e, "%.12e\n", value)
		end
	end
	
	open("n_f.dat", "w") do f_n_f
		for row in eachrow(n_f)
			@printf(f_n_f, "%.12e\t%.12e\t%.12e\n", row...)
		end
	end
	
	open("centroid_faces.dat", "w") do f_centroid_faces
		for row in eachrow(centroid_faces)
			@printf(f_centroid_faces, "%.12e\t%.12e\t%.12e\n", row...)
		end
	end
	
	open("centroid_edges.dat", "w") do f_centroid_edges
		for row in eachrow(centroid_edges)
			@printf(f_centroid_edges, "%.12e\t%.12e\t%.12e\n", row...)
		end
	end
	
	open("n_f_e.dat", "w") do f_n_f_e
		for row in eachrow(n_f_e)
			@printf(f_n_f_e, "%.12e\t%.12e\t%.12e\n", row...)
		end
	end
	
	open("n_fp_e.dat", "w") do f_n_fp_e
		for row in eachrow(n_fp_e)
			@printf(f_n_fp_e, "%.12e\t%.12e\t%.12e\n", row...)
		end
	end
	
	
end

using LinearAlgebra, Printf, DelimitedFiles

function create_edge_coordinates(filename_in::String, filename_out::String, folder_path::String, M::Float64)
    # This function generates files to represent the edges of a polyhedron.
    # The first two components show the vertices forming the edge,
    # and the last two components show the associated faces.
    
    # Define input and output file paths
    file_in = joinpath(folder_path, filename_in)
    file_out = joinpath(folder_path, filename_out)
    
    data = readdlm(file_in, '\n', String)  # Read lines as strings
    
    # Parse the file manually to handle multiple spaces
    vertices = []
    faces = []
    
    for line in data
        tokens = split(line, r"\s+")  # Split on arbitrary whitespace using a regex
        if tokens[1] == "v"
            push!(vertices, parse.(Float64, tokens[2:end]))
        elseif tokens[1] == "f"
            push!(faces, parse.(Int, tokens[2:end]))
        end
    end
    
    # Convert to matrices
    vertices = collect(hcat(vertices...)')
    faces = collect(hcat(faces...)')
    
    # Function to adjust polyhedron data
    vertices = polyhedron_correction(vertices, faces, M, filename_in)
    
    # Determine the number of edges (using an upper bound for non-convex polyhedra)
    qtd_edges = (size(vertices, 1) + size(faces, 1)) * 2
    e_alfa = zeros(Int, qtd_edges, 4)
    
    # Progress bar setup
    println("Finding the coordinates of the edges of the polyhedron...")
    
    s_antes = 0
    for i in 1:size(faces, 1) # Iterate through each face
        AA = Int[]
        nff = count(x -> x[4] == i, eachrow(e_alfa)) # Contact faces for the current face
        
        for j in i:size(faces, 1) # From the current face onward
            int = intersect(faces[i, :], faces[j, :]) # Find contact faces
            if length(int) == 2
                push!(AA, [int..., i, j])
            end
            
            # Exit loop if 3 contact faces are found
            if length(AA) + nff == 3
                break
            end
        end
        
        # Update rows in e_alfa
        s_depois = length(AA) + s_antes
        if !isempty(AA)
            e_alfa[s_antes + 1:s_depois, :] .= AA
        end
        s_antes = s_depois
    end
    
    # Remove rows with all zeros (excess rows)
    e_alfa = e_alfa[any(e_alfa .!= 0, dims=2), :]
    
    # Write the edge data to the output file
    open(file_out, "w") do file
        for row in eachrow(e_alfa)
            @printf(file, "%d %d %d %d\n", row...)
        end
    end
    
    # Write vertices and faces to separate files
    open(joinpath(folder_path, "faces.dat"), "w") do file
        for row in eachrow(faces)
            @printf(file, "%.20e %.20e %.20e\n", row...)
        end
    end
    
    open(joinpath(folder_path, "vertex.dat"), "w") do file
        for row in eachrow(vertices)
            @printf(file, "%.20e %.20e %.20e\n", row...)
        end
    end
end


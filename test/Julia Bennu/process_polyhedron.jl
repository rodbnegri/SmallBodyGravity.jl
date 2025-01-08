include("polyhedron_correction.jl")
include("create_edge_coordinates.jl")

function process_polyhedron(filename::String)
    # Ensure the file exists in the current working directory
    if !isfile(filename)
        throw(ArgumentError("File $filename does not exist in the current directory or specified path."))
    end
    
    # Use the current working directory as folder_path
    folder_path = pwd()
    
    # Call the create_edge_coordinates function
    create_edge_coordinates("bennu.obj", "Edge_Coords.dat", folder_path, 7.329e10)
end

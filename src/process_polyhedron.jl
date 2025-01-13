

"""
    process_polyhedron(filename::String, Mass::Float64, folder_path::String = pwd())

Processes a polyhedron file by generating edge coordinates and additional data for gravity modeling.

# Arguments
- `filename::String`: Name of the polyhedron file to process. The file must be in the Wavefront format with only the vertex ("v"s) and faces ("f"s) in the file
- `Mass::Float64`: Mass of the small body in kilograms.
- `folder_path::String`: Path to the folder containing the input file. Defaults to the current working directory (`pwd()`).

"""

function process_polyhedron(filename::String, Mass::Float64, folder_path::String = pwd()) 
    
    # Construct the full file path 
    file_path = joinpath(folder_path, filename)
    # Ensure the specified file exists
    if !isfile(file_path)
        throw(ArgumentError("File $filename does not exist in the specified path: $folder_path."))
    end
    
    # Step 1: Generate edge coordinates
    # This function takes the input file, processes it, and generates "Edge_Coords.dat" in the specified folder.
    create_edge_coordinates(file_path, "Edge_Coords.dat", folder_path, Mass)
    
    # Step 2: Further process the polyhedron and edge data
    # This function uses the input file and the generated edge data to compute additional information.
    polyhedron_data(file_path, "Edge_Coords.dat")
end

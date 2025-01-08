using LinearAlgebra, Printf

function polyhedron_correction(vertex_old::Matrix{Float64}, faces::Matrix{Int}, M::Float64, filename_in::String)
    # Models from: Dobrovolskis - Inertia of Any Polyhedron
    
    # Finding volume and corresponding density
    Volume = 0.0
    for simplice in 1:size(faces, 1)
        vertex_1 = vertex_old[faces[simplice, 1], :]
        vertex_2 = vertex_old[faces[simplice, 2], :]
        vertex_3 = vertex_old[faces[simplice, 3], :]
        Volume += dot(cross(vertex_2, vertex_3), vertex_1) / 6.0
    end
    sigma = M / Volume
   
  # Finding the center of mass
    R = zeros(3)
    for simplice in 1:size(faces, 1)
        vertex_1 = vertex_old[faces[simplice, 1], :]
        vertex_2 = vertex_old[faces[simplice, 2], :]
        vertex_3 = vertex_old[faces[simplice, 3], :]
        dV = dot(cross(vertex_2, vertex_3), vertex_1) / 6.0
        dR = (vertex_1 + vertex_2 + vertex_3) / 4.0
        R += dV * dR / Volume
    end
    
    # Finding the products of inertia and the inertia matrix
    P = zeros(3, 3)
    for simplice in 1:size(faces, 1)
        for j in 1:3
            for k in 1:3
                D = vertex_old[faces[simplice, 1], :]
                E = vertex_old[faces[simplice, 2], :]
                F = vertex_old[faces[simplice, 3], :]
                dV = dot(cross(E, F), D) / 6.0
                P[j, k] += sigma * dV / 20.0 * (
                2 * D[j] * D[k] + 2 * E[j] * E[k] + 2 * F[j] * F[k] +
                D[j] * E[k] + D[k] * E[j] + D[j] * F[k] + D[k] * F[j] +
                E[j] * F[k] + E[k] * F[j]
                )
            end
        end
    end
    
    I = [P[2, 2] + P[3, 3] -P[1, 2] -P[1, 3];
    -P[1, 2] P[1, 1] + P[3, 3] -P[2, 3];
    -P[1, 3] -P[2, 3] P[1, 1] + P[2, 2]]
    
    # Parallel Axis Theorem: find inertia tensor relative to center of mass
    ICM = I - M * [R[2]^2+R[3]^2 -R[1]*R[2] -R[1]*R[3];
    -R[1]*R[2] R[1]^2+R[3]^2 -R[2]*R[3];
    -R[1]*R[3] -R[2]*R[3] R[1]^2+R[2]^2] 
    
    # Eigenvalues and eigenvectors
    eig_val, eig_vec = eigen(ICM)
    eig_vec[diagm(eig_val) .< 0] .= -eig_vec[diagm(eig_val) .< 0] # eigen vector matrix diagonal should be positive to form the new coordinate system
    
    # Aligning with the principal axes
    vertex = zeros(size(vertex_old))  
    for simplice in 1:size(vertex_old, 1)
        vertex[simplice, :] = (vertex_old[simplice, :] .- R)' * eig_vec
    end
    
    # Equivalent ellipsoid
    A, B, C = eig_val[1], eig_val[2], eig_val[3]
    a = sqrt(5 * (B + C - A) / 2 / M)
    b = sqrt(5 * (A + C - B) / 2 / M)
    c = sqrt(5 * (A + B - C) / 2 / M)
    
    
    # Finding new properties
    
    # Initialize the center of mass
    Rnew = zeros(3)
    
    # Loop through each simplex (face)
    for simplice in 1:size(faces, 1)
        vertex_1 = vertex[faces[simplice, 1], :]
        vertex_2 = vertex[faces[simplice, 2], :]
        vertex_3 = vertex[faces[simplice, 3], :]
        
        dV = dot(cross(vertex_2, vertex_3), vertex_1) / 6  # Volume of the current simplex
        dR = (vertex_2 + vertex_3 + vertex_1) / 4         # Center of mass of the current simplex
        Rnew = dV .* dR ./ Volume + Rnew
    end
    
    # Initialize the inertia products matrix
    P = zeros(3, 3)
    
    # Loop to compute the products of inertia
    for simplice in 1:size(faces, 1)
        for j in 1:3
            for k in 1:3
                D = vertex[faces[simplice, 1], :]
                E = vertex[faces[simplice, 2], :]
                F = vertex[faces[simplice, 3], :]
                
                dV = dot(cross(E, F), D) / 6  # Volume of the current simplex
                
                P[j, k] += sigma * dV / 20 * (
                2 * D[j] * D[k] + 2 * E[j] * E[k] + 2 * F[j] * F[k] +
                D[j] * E[k] + D[k] * E[j] +
                D[j] * F[k] + D[k] * F[j] +
                E[j] * F[k] + E[k] * F[j]
                )
            end
        end
    end
    
    # Compute the inertia matrix
    Inew = [P[2, 2]+P[3, 3] -P[1, 2] -P[1, 3];
    -P[1, 2] P[1, 1]+P[3, 3] -P[2, 3];
    -P[1, 3] -P[2, 3] P[1, 1]+P[2, 2]]
    
    
    # Write results to a file
    open("polyhedron_properties.txt", "w") do file
        write(file, "Number of faces: $(size(faces, 1)) \n\n")
        write(file, "Mass [kg]: $(M) \n\n")
        write(file, "Volume [km^3]: $(Volume) \n\n")
        write(file, "Density [kg/m^3]: $(sigma * 1e-9) \n\n")
        # Write inertia matrix row by row
        write(file, "Inertia matrix [kg.km^2]:\n")
        for row in 1:size(Inew, 1)
            write(file, "$(Inew[row, :])\n")  # Write each row on a new line
        end
        write(file, "\n")
        write(file, "Center of mass [km]: $(Rnew)\n\n")
        write(file, "Equivalent ellipsoid:\n")
        write(file, "a x b x c [km]: $(a) x $(b) x $(c)\n")
        write(file, "Equivalent volume [km^3]: $(4 / 3 * π * a * b * c)\n")
        write(file, "Equivalent density [kg/m^3]: $(3 * M / (4 * π * a * b * c) * 1e-9)\n\n")
        write(file, "***OLD PROPERTIES*** as found in $(filename_in):\n")
        write(file, "Inertia matrix [kg.km^2]:\n")
        for row in 1:size(I, 1)
            write(file, "$(I[row, :])\n")  # Write each row on a new line
        end
        write(file, "\n")
        write(file, "Center of mass [km]: $(R)\n")
    end
    
    # Remove the extension from the filename_in
    filename_base = splitext(filename_in)[1]
    
    # Construct the new filename
    filename_out = "$(filename_base)_corrected.obj"
    
    # Make new .obj file
    open(filename_out, "w") do file
        # Write vertex coordinates, each row as a 'v' line
        for i in 1:size(vertex, 1)
            write(file, "v\t$(vertex[i, 1])\t$(vertex[i, 2])\t$(vertex[i, 3])\n")
        end
        
        # Write faces, each row as an 'f' line
        for j in 1:size(faces, 1)
            write(file, "f\t$(faces[j, 1])\t$(faces[j, 2])\t$(faces[j, 3])\n")
        end
    end
    
    
    return vertex
end


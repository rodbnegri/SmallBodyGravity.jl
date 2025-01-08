module SmallBodyGravity

using Pkg
Pkg.add("Printf")
Pkg.add("LinearAlgebra")
Pkg.add("DelimitedFiles")

using Printf
using LinearAlgebra, Printf, DelimitedFiles

include("spherical_harmonics.jl")
include("polyhedron_model.jl")
include("polyhedron_correction.jl")
include("create_edge_coordinates.jl")

export spherical_harmonics 
export polyhedron_model
export polyhedron_correction
export create_edge_coordinates

end


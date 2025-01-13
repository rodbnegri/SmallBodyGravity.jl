
# SmallBodyGravity


[![Build Status](https://github.com/rodbnegri/SmallBodyGravity.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/rodbnegri/SmallBodyGravity.jl/actions/workflows/CI.yml?query=branch%3Amaster) ![Status](https://img.shields.io/badge/status-WIP-yellow)

⚠️ **This package is in development. Some features may not be fully functional, and parts of the code have not been tested.**

This Julia package provides the equations to model the gravity field of small bodies. The current implementation focuses on the mathematical equations for modeling the gravity field, that I originally used in the papers:

[1] **Negri, R.B. and Prado, A.F., 2022. Autonomous and Robust Orbit-Keeping for Small-Body Missions. Journal of Guidance, Control, and Dynamics, 45(3), pp.587-598. DOI: 10.2514/1.G005863**

[2] **Negri, R.B., Prado, A.F., Chagas, R.A. and Moraes, R.V., 2024. Autonomous Rapid Exploration in Close-Proximity of Asteroids. Journal of Guidance, Control, and Dynamics, 47(5), pp.914-933. DOI: 10.2514/1.G007186**

, and in my Ph.D. Thesis:

[3] **Negri, R. B. (2022). Estudo de dinâmica, guiagem, navegação e controle aplicado à deflexão de asteroides. [Doctoral thesis, Instituto Nacional de Pesquisas Espaciais (INPE)]. Instituto Nacional de Pesquisas Espaciais. http://urlib.net/ibi/8JMKD3MGP3W34T/46L9G6P**

If you use this code for academic purposes, please consider citing one of the references above. This would encourage me to share more of my work. For example: "We used the SmallBodyGravity package in Julia, originally developed by Dr. Negri for the works [1, 2, 3]."

## Package Overview

The `process_polyhedron` function pre-processes a polyhedron file in the Wavefront .obj format, which should only contain "v" (vertices) and "f" (faces) lines (other types of data are not currently supported). The function aligns the coordinates of the polyhedron with the principal axes of the small body, based on the specified mass, and generates the necessary data for gravity field calculations. 

Additionally, the function generates a `polyhedron_properties.txt` file, which includes properties of the small body (e.g., volume, density, equivalent ellipsoid, etc.). A corrected `.obj` file is also produced.

### After Running `process_polyhedron`

Once the `process_polyhedron()` function has been executed (you just need to execute it a single time), the generated `.dat` files (such as `centroid_edges.dat`, `e_e.dat`, etc.) can be used to calculate the gravity field. To do this, include each of the `.dat` files using the `readdlm()` function, and add them to the parameter tuple `p` in `polyhedron_model()`.

Refer to the `example.jl` file located in the `example` folder for a full example of usage.



## Future Work

A MATLAB code that calculates the harmonic coefficients from the small body's polyhedron is being transitioned to Julia and it should be available in the coming months.

## Installation

To install the package, you can use the following command in Julia:

```julia
using Pkg
Pkg.add("SmallBodyGravity")

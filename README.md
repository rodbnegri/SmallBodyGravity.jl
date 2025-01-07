
# SmallBodyGravity

[![Build Status](https://github.com/rodbnegri/SmallBodyGravity.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/rodbnegri/SmallBodyGravity.jl/actions/workflows/CI.yml?query=branch%3Amaster)

This Julia package provides the equations to model the gravity field of small bodies. The current implementation focuses on the mathematical equations for modeling the gravity field, as described in the papers:

- **Batista Negri, R. and Prado, A.F., 2022. Autonomous and Robust Orbit-Keeping for Small-Body Missions. Journal of Guidance, Control, and Dynamics, 45(3), pp.587-598.**
- **Negri, R.B., Prado, A.F., Chagas, R.A. and Moraes, R.V., 2024. Autonomous Rapid Exploration in Close-Proximity of Asteroids. Journal of Guidance, Control, and Dynamics, 47(5), pp.914-933.**

## Package Overview

At this stage, the package includes the core equations for the gravity field model of small bodies. The polyhedron model of the small body is represented, but please note that the coordinates of the polyhedron must be aligned as follows:
- The origin should be at the center of mass of the small body.
- The axes should be aligned with the principal axes of inertia.

## Future Work

Currently, the code to calculate the center of mass and align the axes of inertia is implemented in MATLAB. I am in the process of transitioning this MATLAB code to Julia, and this functionality will be integrated into the package in the future.

## Installation

To install the package, you can use the following command in Julia:

```julia
using Pkg
Pkg.add("SmallBodyGravity")
.yml/badge.svg?branch=master)](https://github.com/rodbnegri/SmallBodyGravity.jl/actions/workflows/CI.yml?query=branch%3Amaster)

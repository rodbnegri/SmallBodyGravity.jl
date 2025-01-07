
# SmallBodyGravity


[![Build Status](https://github.com/rodbnegri/SmallBodyGravity.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/rodbnegri/SmallBodyGravity.jl/actions/workflows/CI.yml?query=branch%3Amaster) ![Status](https://img.shields.io/badge/status-WIP-yellow)

⚠️ **This package is in development. Some features may not be fully functional, and parts of the code have not been tested.**

This Julia package provides the equations to model the gravity field of small bodies. The current implementation focuses on the mathematical equations for modeling the gravity field, that I originally used in the papers:

- **Negri, R.B. and Prado, A.F., 2022. Autonomous and Robust Orbit-Keeping for Small-Body Missions. Journal of Guidance, Control, and Dynamics, 45(3), pp.587-598. DOI: 10.2514/1.G005863**
- **Negri, R.B., Prado, A.F., Chagas, R.A. and Moraes, R.V., 2024. Autonomous Rapid Exploration in Close-Proximity of Asteroids. Journal of Guidance, Control, and Dynamics, 47(5), pp.914-933. DOI: 10.2514/1.G007186**

## Package Overview

At this stage, the package includes the core equations for the gravity field model of small bodies. The polyhedron model of the small body is represented, but please note that the coordinates of the polyhedron must be aligned as follows:
- The origin should be at the center of mass of the small body for the chosen mass.
- The axes should be aligned with the principal axes of inertia for the chosen mass.

## Future Work

Currently, the code to calculate the center of mass and align the axes of inertia for the chosen mass is implemented in MATLAB. I am in the process of transitioning this MATLAB code to Julia, and this functionality will be integrated into the package in the future.

A MATLAB code that calculates the multipolar expansion coefficients from the polyhedron is also being transitioned to Julia.

## Installation

To install the package, you can use the following command in Julia:

```julia
using Pkg
Pkg.add("SmallBodyGravity")

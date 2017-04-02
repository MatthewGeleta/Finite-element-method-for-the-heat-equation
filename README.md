# Finite-element-method-for-the-heat-equation
This repository contains MATLAB code for a finite element solution to the stochastic heat equation with non-zero
Dirichlet boundary conditions and forcing function on a non-simple domain.

The code is based on the paper "Semigroups and Finite Elements for the Stochastic Heat Equation", by Matthew Geleta,
submitted as a special topic for an MSc degree in Mathematics at The University of Oxford.

Mathematical details:
Galerkin approximation using triangular cells and Lagrange basis elements.
Time-stepping is implemented using a backward-Euler scheme.
Stochastic integrals are computed using the Euler-Maruyama method.

Files:
FEM_for_heat_equation - deterministic case.
FEM_for_the_Stochastic_Heat_Equation - stochastic case.
Generate_FEM_Mesh - triangulation of a non-simple domain.
Dloc - function for local mass matrix.
Bloc - function for local stiffness matrix.
FEM_RHS_t - forcing function for heat equation.
Other files are dependencies for Generate_FEM_Mesh.

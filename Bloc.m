function B = Bloc(coords)
%{
Bloc returns the local stiffness matrix B determined by integrating
derivatives of Lagrange basis elements over a triangle with coordinate
vertices given by the input coords.

Input:
coords:: 3-by-2 matrix of vertex coordinates.
    First column are x-coords, second are y-coords.

Output:
D:: local matrix

Author: Matthew Geleta
Date: 02/04/2017

%}

G = [ones(1,3); coords']\[0,0;1,0;0,1];
B = det([ones(1,3);coords'])*G*G'/2;
end
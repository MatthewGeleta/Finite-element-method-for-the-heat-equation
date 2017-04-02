function D = Dloc(coords)
%{
Dloc computes the local contribution to the global mass matrix of the 
triangular cell with vertex coordinates 'coords', using Lagrange
basis elements.

Input:
coords:: 3-by-2 matrix of vertex coordinates.
    First column are x-coords, second are y-coords.
Output:
D:: local matrix

Author: Matthew Geleta
Date: 02/04/2017
%}
D = det([1,1,1; coords'])*[2,1,1;1,2,1;1,1,2]/24;
end
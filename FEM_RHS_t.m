function b = FEM_RHS_t(p,t,n,dt)

M = size(t,1); % Number of triangles
N = size(p,1);
b = zeros(N,1);

% Comment out the for loop for case of zero external forcing:
%{
% Loop over elements
for j = 1:M 
    nodes = t(j,:);
    coords = p(nodes,:);
    b(nodes) = b(nodes) + ...
            det([1,1,1; p(t(j,:),:)']) * ...
            dt*f_forcing( sum(p(t(j,:), :))/3 )/6; % approximate the integral using
                % Simpson's rule with Barycentric centre evaluation for f
end
%}

end
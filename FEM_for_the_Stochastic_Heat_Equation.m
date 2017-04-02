%% Stochasitc heat equation by finite elements
%{
Author: Matthew Geleta
Date: 02/04/2017

Overview:
The following script simulates the evolution of a two-dimensional 
stochastic heat equation in a non-simple domain using the Galerkin finite 
element method for spatial the discretisation, and an Euler-Maruyama method
for stochastic integration.
%}

clear all;

%% Offline computations
%
%% Mesh generation
%
Generate_FEM_Mesh;
%}
%% Set parameters
%
dt = 0.01; % time-step
sd = sqrt(dt); % standard deviation of the Brownian noise increments
scale = 0.01; % magnitude of noise term
Tstep = 20; % maximum time-steps
Tmax = Tstep*dt; % maximum simulation time
N = size(p,1); % number of mesh nodes
M = size(t,1); % number of cells in triangulation
%}
%% Get boundary and interior nodes
%
AllNodes = 1:N; % list of all nodes

% Nodes for each inner circular boundary:
B1 = (sum( (p-[p1x,p1y]).^2, 2) < R1^2+h0/100 );
B2 = (sum( (p-[p2x,p2y]).^2, 2) < R2^2+h0/100 );
B3 = (sum( (p-[p3x,p3y]).^2, 2) < R3^2+h0/100 );

% Nodes for large circular boundary:
Bout = (sum( p.^2, 2) > Rout^2-h0/100 );

HotNodes = (B1|B2|B3); % nodes on hot (inner) boundaries
ColdNodes = (Bout); % nodes on cold (outer) boundary
Bdry = AllNodes(HotNodes | ColdNodes); % all boundary nodes
IN = setdiff(AllNodes, Bdry); % all interior nodes
%}
%% Visualise boundary nodes
%
figure()
scatter(p(B1,1),p(B1,2))
hold on
scatter(p(B2,1),p(B2,2))
scatter(p(B3,1),p(B3,2))
scatter(p(Bout,1),p(Bout,2))
title('Boundary nodes')
legend('Boundary 1', 'Boundary 2', 'Boundary 3', 'Outer boundary')
%}
%% Random boundary temperatures
%
% Generate boundary temperatures as uniformly random on [0.5,1.5]:
T1 = 0.5 + rand; % temperature for B1 boundary
T2 = 0.5 + rand; % temperature for B2 boundary
T3 = 0.5 + rand; % temperature for B3 boundary
Tout = 0; % temperature for outer boundary
%}
%% Assembly of global mass and stiffness matrices
%
D = sparse(N,N); % global mass matrix
B = sparse(N,N); % global stiffness matrix

% Loop over cells in mesh, compute local contributions using exact
% quadrature rules, and add to global matrices:
for j = 1:M
   nodes = t(j,:); % vertices of triangles
   coords = p(nodes,:); % coordinates of vertices
   % Compute D by quadrature rule of Lagrange basis
   D(nodes,nodes) = D(nodes,nodes) + Dloc(coords);
   % Compute B by quadrature rule for derivatives of Lagrange basis
   B(nodes,nodes) = B(nodes,nodes) + Bloc(coords);
end
%}
%% Array initialisation
%
c = zeros(N,Tstep); % coefficient vector at each time-step
% Inner boundaries:
c(B1) = T1;
c(B2) = T2;
c(B3) = T3;
% Outer boundary:
c(ColdNodes,1) = Tout;
%}

%% Online computations
%
%% Commence time-stepping
%
for n = 1:Tstep-1
    % Temperature-dependent stochastic noise term for this time-step
    dw = sd*randn(N,1).*c(:,n)*scale;
    % Forcing function at time n*dt
    b = FEM_RHS_t(p,t,n,dt);
    % Contribution from previous time-step:
    b = b + D * c(:,n);
    % Dirichlet boundary conditions:
    cn = sparse(N,1);
    cn(B1) = T1;
    cn(B2) = T2;
    cn(B3) = T3;
    cn(ColdNodes) = Tout;
    % Modify RHS using boundary conditions:
    b = b - (dt * B + D) * cn;
    % Incorparate noise via Euler-Maruyama step:
    b = b + dw;
    % Time-step (only update interior nodes):
    cn(IN) = (D(IN,IN) + dt*B(IN,IN))\b(IN);
    % Store solution for current time-step
    c(:,n+1) = cn;
end
%}

%% Plot individual time-step
%
fg1 = figure(1)
set(fg1, 'Position', [5 5 1000 1000])
n = 5; % time-step to plot
trisurf(t,p(:,1),p(:,2),c(:,n),'facecolor','interp')
%pbaspect([1 1 1])
%set(gca,'visible','off');
hold off;
title('Heat distribution')
zlabel('Temperature')
%}
%% Plot evolution of heat flow
%
figure()
for n = 1:Tstep
trisurf(t,p(:,1),p(:,2),c(:,n),'facecolor','interp')
zlim([-0.2,1.2]);
pause(0.1)
end
%}



%% Run multiple stochastic simulations and compute expectation
%
sims = 10; % number of simulations
coefs = zeros(N,sims);

for m = 1:sims
    
    T1 = 0.5 + rand;%0.75; % temperature for B1 boudaries
    T2 = 0.5 + rand;%0.5; % temperature for B2 boundaries
    T3 = 0.5 + rand;%1; % temperature for B3 boundaries
    
    c = zeros(N,Tstep); % coefficient vector at each time-step
    %c(HotNodes,1) = 1; % hot boundary
    c(B1) = T1;
    c(B2) = T2;
    c(B3) = T3;
    c(ColdNodes,1) = 0; % cold boundary
    
    for n = 1:Tstep-1
    % Noise term for this time step
    dw = sd*randn(N,1).*c(:,n)*scale;
    % Forcing function at time n*dt
    b = FEM_RHS_t(p,t,n,dt);
    % Contribution from previous time-step:
    b = b + D * c(:,n);
    % Dirichlet boundary conditions:
    cn = sparse(N,1);
    %cn(HotNodes) = 1;
    cn(B1) = T1;
    cn(B2) = T2;
    cn(B3) = T3;
    cn(ColdNodes) = 0;
    % Modify RHS using boundary conditions:
    b = b - (dt * B + D) * cn;
    % Modify RHS using Euler-Marayuama noise term:
    b = b + dw;
    % Time-step using Euler step (only consider interior nodes):
    cn(IN) = (D(IN,IN) + dt*B(IN,IN))\b(IN);
    % Store solution for current time-step
    c(:,n+1) = cn;
end
    coefs(:,m) = c(:,end);
    
end

% Compute expectation
C = mean(coefs,2);

fg1 = figure(1)
set(fg1, 'Position', [5 5 1000 1000])
trisurf(t,p(:,1),p(:,2),C,'facecolor','interp')
%pbaspect([1 1 1])
 set(gca,'visible','off');
 hold off;
%title('Heat distribution: linear Lagrange basis')
%zlabel('Temperature')




%
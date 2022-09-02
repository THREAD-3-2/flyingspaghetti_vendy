%% Input data Free flight of a beam

n_ele = 2;                                      % number of elements

index = 1:n_ele+1;                              % indices for the nodes
xc = 0:6/n_ele:6;                               % initial configuration
zc = 8:-8/n_ele:0;

nodal_c = [index' xc' zeros(n_ele+1,1) zc'];    %  nodal coordinates

% number of nodes
n_node = size(nodal_c,1);

% nodal connections
pq = [index(1:end-1)' index(1:end-1)' index(2:end)'];

%number of elements
n_elem =  size(pq,1);

%% Material data

% modulus of elasticity
E = 10^4;

%shear modulus
G = 10^4;          

% Shear area A1 [m2]:
A1 = 1;  

% Shear area A2 [m2]:
A2 = 1; 

% Area of cross section A [m2]:
A = 1; 

% Moment of inertia about 1 [m4]:
J1 = 0.05; 

% Moment of inertia about 2 [m4]:
J2 = 0.05;   

% Torsional moment of cross-section [m4]:
Jt = 0.05;

% Density
rho = 1;

% Inertial moment
Jp = diag([10 10 10]);

% Cross-sectional tangent modulus
C11 = diag([1 E*A G*A1 G*A2]);
C12 = zeros(4);
C21 = zeros(4);
C22 = diag([1 G*Jt E*J1 E*J2]);

%no initial strains and the beam cross-sections are initially orthogonal to
%the axis of the beam
gama0 = [-1;0;0];

%% Order of the element
n_order = 3;
ng = n_order+1;             % number of integration points
[Xg,Wg] = GaussInt(ng);     % Gaussian points and weights for numerical integration
if n_order == 1
    N_node = size(nodal_c,1);
else
    N_node = n_order*n_elem+1;
end

%Convergence criteria for the termination of Newton iteration
convergence_criteria = 10^-8;

% Nodal coordinates for higher order elements
INDEX = 1:1:n_elem*n_order+1;
XC = 0:6/(n_order*n_ele):6;
ZC = 8:-8/(n_ele*n_order):0;
Nodal_C = [INDEX' XC' zeros(n_elem*n_order+1,1) ZC'];  

%% time intervals
dt = 0.1;               % time interval
total_time = 30;        % total time
t_i = 0;                % time initialization
%% Point loads and moments


% F_ele = [nodal_c(i) fx fy fz Mx My Mz;
%          nodal_c(j) fx fy fz Mx My Mz;
%           .....]
F_ele = [n_node 1/10 0 0 0 1 1/2];
% F_node = [Nodal_C(i) fx fy fz Mx My Mz;
%           Nodal_C(j) fx fy fz Mx My Mz;
%           .....]
F_node = [N_node 1/10 0 0 0 1 1/2];

% Boundary condition
% Bc = [i fix(1) fix(2)..... fix(6);
%       j fix(1) fix(2)..... fix(6)
%       .....]
Bc = [];                                        % Boundary condition for a free beam

%% Input data Imbrahimbeg_test

n_ele = 50;              % number of elements

Len = 10;                % Length of the elements
index = 1:n_ele+1;                % indices for the nodes
interval = Len*(index-1)/n_ele;              % Intervals

nodal_c = [index' interval' zeros(n_ele+1,2)];    %  nodal coordinates

% number of nodes
n_node = size(nodal_c,1);


%% element connection

i = 1:n_ele;

% nodal connections
pq = [i' i' (i+1)'];
%number of elements
n_elem =  size(pq,1);

ele_type = 0;                               %(1 - curved, 0 - straight)

strait_ele = repmat([1 0 1 0 0], n_ele+1,1);

% Bc = [i fix(1) fix(2)..... fix(6)]


Bc = [1 1 2 0 4 5 0];      
% Bc = [];


relaxation = {};


%% Material data
%material properties
material = 1;           % linear elastic model

% modulus of elasticity
E = 10^4;

% poissons ratio
%nu = 0;

%shear modulus
G = E;                      %E/(2*(1+nu));

% Shear area A1 [m2]:
A1 = 1;  

% Shear area A2 [m2]:
A2 = 1; 

% Area of cross section A [m2]:
A = 1; 

% Moment of inertia about 1 [m4]:
J1 = 0.1; 

% Moment of inertia about 1 [m4]:
J2 = 0.1;   

% Torsional moment of cross-section [m4]:
Jt = 0.1;

rho = 1;

% Inertial moment

Jp = diag([20 10 10]);

C11 = diag([1 E*A G*A1 G*A2]);
C12 = zeros(4);
C21 = zeros(4);	
C22 = diag([1 G*Jt E*J1 E*J2]);

gama0 = [-1;0;0];

%% integration points
n_order = 1;
ng = n_order+1;
ngr = n_order+1;

int = 'G';
if int == 'G'

    [Xgr,Wgr] = GaussInt(ngr);
    [Xg,Wg] = GaussInt(ng);

elseif int == 'L'
    
    [Xgr,Wgr] = LobattoInt(ngr);
    [Xg,Wg] = LobattoInt(ng);

end

if n_order == 1
    
    N_node = size(nodal_c,1);
    
else
    
    N_node = n_order*n_elem+1;
    
end
    

%Newton iteration
convergence_criteria = 10^-8;

%Max iteration
max_it = 18;

INDEX = 1:1:n_elem*n_order+1;
INTERVAL = Len*(INDEX-1)/(n_ele*n_order);
Nodal_C = [INDEX' INTERVAL' zeros(n_elem*n_order+1,2)];  

%% time intervals
dt = 0.1;
total_time = 10;
t_i = 0;

%% Load
             %node 
%nodal_mass = [i    M   Tx  Ty  Tz]


nodal_mass = [];
Load = 80;
F_ele = [1 0 0 0.05*Load 0 0 Load];
F_node = [1 0 0 0.05*Load 0 0 Load];


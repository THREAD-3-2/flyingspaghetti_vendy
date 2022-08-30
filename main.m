%% main
clc;
clear all;
close all;
                                    % #########################################################
                                    % ##         VELOCITY BASED FINITE ELEMENT FORMULATION   ##
                                    % ##                                                     ##
                                    % ##                                                     ##
                                    % ##            SEE THE PAPER FOR MORE DETAILS           ##
                                    % #########################################################


% [1] Zupan, E. and Zupan, D., 2018. On conservation of energy and kinematic compatibility in dynamics of nonlinear velocity-based three-dimensional beams....
% Nonlinear Dynamics, 95(2), pp.1379-1394.
fprintf('\n Velocity based elements \n');

%% input data

Flying_spaghetti                    % Input data for the free flight of a flexible beam

%% Preprocessing
                                    % ***************************************************
                                    % *  READING THE INPUT AND PREPROCESSING THE DATA   *
                                    % ***************************************************

fprintf('\n * Number of nodes: %i \n', n_node);
fprintf('\n * Number of elements: %i \n', n_elem);
fprintf('\n * Time step size: %i \n', dt);
% Supports
fix = [];

for i=1:size(Bc,1)
    fix = union(fix,nonzeros(Bc(i,2:size(Bc,2)))+6*(Bc(i,1)-1));
end

dof_free = 1:6*N_node;              % number of degrees of freedom
dof_free(fix) = [];                 % free degrees of freedom 


%% Length of elements
elem_l = 0;
elem_len = zeros(n_node,1);
ele_c = zeros(n_elem,n_order+1);
for i = 1:n_elem
    elem_l = sqrt((nodal_c(i+1,2) - nodal_c(i,2))^2 + (nodal_c(i+1,3) - nodal_c(i,3))^2 + (nodal_c(i+1,4) - nodal_c(i,4))^2);
    elem_len(i+1,1) = elem_len(i,1) + elem_l; 
    vec_el(1:3,i) = [nodal_c(pq(i,3),2) - nodal_c(pq(i,2),2), nodal_c(pq(i,3),3) - nodal_c(pq(i,2),3), nodal_c(pq(i,3),4) - nodal_c(pq(i,2),4)]';
    for j = 1:n_order+1
        if n_order == 1
            ele_c(i,j) = elem_len(i)+elem_len(j);
        else
            ele_c(i,j) = elem_len(i) + (j-1)*((elem_len(i+1)-elem_len(i))/(n_order));
        end
    end
end

%% Initial rotations

q_identity = [1;0;0;0];
rot_ini = [0.6,0,0.8;0,1,0;-0.8,0,0.6];
thta_ini = Spurrier(rot_ini); 

%% Initialization at time 0

time_n = 1;
T_energy(1,time_n) = 0;
pt(time_n) = 0;
v{time_n} = zeros(3,N_node);
omega{time_n} = zeros(3,N_node);
v_bar{time_n} = zeros(3,N_node);
omega_bar{time_n} = zeros(3,N_node);
gama{time_n} = zeros(4,N_node);
kapa{time_n} = zeros(4,N_node);
rg{time_n} = zeros(3,N_node);
rgdx{time_n} = zeros(3,N_node);
qn{time_n} = zeros(4,N_node);
kn{time_n} = zeros(4,N_node);

knI{time_n} = zeros(4,ng,n_elem);
rgdx_I{time_n} = zeros(3,ng,n_elem);
gamaI{time_n} = zeros(4,ng,n_elem);
kapaI{time_n} = zeros(4,ng,n_elem);
vI{time_n} = zeros(3,ng,n_elem);
omegaI{time_n} = zeros(3,ng,n_elem);

for e = 1:n_elem
    q0_R(:,e) = quaternion_rep(thta_ini);
    ind_eL = n_order*(e-1) + 1;
    ind_eR = n_order*e + 1;   
    for ni = ind_eL:ind_eR
        kn{1}(:,ni) = q_identity;
        qn{1}(:,ni) = Quaternion_product(q0_R(:,e),q_identity);
        rgdx{1}(:,ni) = -rot_ini*gama0;
        rg{1}(:,ni) = Nodal_C(ni,2:4);  
    end
    for ij = 1:ng
        knI{1}(:,ij,e) = q_identity;
        rgdx_I{1}(:,ij,e) = -rot_ini*gama0;
    end    
end

%% Analysis

tic

while t_i<=total_time                                           % Loop over time
    
    % initialization for next time step
    T_energy(1,time_n+1) = 0;
    v{time_n+1} = zeros(3,N_node);
    omega{time_n+1} = zeros(3,N_node);
    v_bar{time_n+1} = zeros(3,N_node);
    omega_bar{time_n+1} = zeros(3,N_node);
    gama{time_n+1} = zeros(4,N_node);
    kapa{time_n+1} = zeros(4,N_node);
    rg{time_n+1} = zeros(3,N_node);
    rgdx{time_n+1} = zeros(3,N_node);
    qn{time_n+1} = zeros(4,N_node);
    
    knI{time_n+1} = zeros(4,ng,n_elem);
    rgdx_I{time_n+1} = zeros(3,ng,n_elem);
    gamaI{time_n+1} = zeros(4,ng,n_elem);
    kapaI{time_n+1} = zeros(4,ng,n_elem);
    vI{time_n+1} = zeros(3,ng,n_elem);
    omegaI{time_n+1} = zeros(3,ng,n_elem);
    
    fprintf('\n time step of iteration: %6i\n',t_i)
    disp(' --------------------------------------------')

    %% Load factor and Load
    FM_elem = zeros(6*(n_elem+1),1);
    FM_node = zeros(6*(N_node),1);
    hy1 = triangular_load(t_i+dt);                                  % triangular load (see example 5.1 in [1] for more details)
    hy0 = triangular_load(t_i);
    hy_bar = (1/2)*(hy1+hy0);
    FM_elem(6*(F_ele(1)-1)+1:6*F_ele(1),1) = dt*hy_bar.*(F_ele(2:7))';         % Point forces and moments in fixed frame
    FM_node(6*(F_node(1)-1)+1:6*F_node(1),1) = dt*hy_bar.*(F_node(2:7))';
    
    fprintf('\n\n ------------------------------------------------------------------------------\n');
    fprintf('%2s %15s %16s\n','iteration','NORM (corrections)','NORM(Residual)');
    disp(' ------------------------------------------------------------------------------')
    
    %% Predictor
    
    v{time_n+1} = v{time_n};
    omega{time_n+1} = omega{time_n};
    v_bar{time_n+1} = (1/2)*(v{time_n+1} + v{time_n});                                      % Predictor step ([V_n+1,Omega_n+1]= [V_n,Omega_n])
    omega_bar{time_n+1} = (1/2)*(omega{time_n+1} + omega{time_n});
        
    vI{time_n+1} = vI{time_n};
    omegaI{time_n+1} = omegaI{time_n};                                               

    %% Newton iteration
    iter = 0;                                                                               % Initiating Newton iterations
    lim_value = 1;                                                                          % norm(residual and corrections) initialized to 1
    while lim_value >= convergence_criteria                                                 % Loop over Newton iteration
    
        [rgdx_I{time_n+1}, vI{time_n+1}, omegaI{time_n+1}, knI{time_n+1}, kapaI{time_n+1}, gamaI{time_n+1}, v_bar{time_n+1},...
        omega_bar{time_n+1}, lim_value, T_energy(time_n+1), f_lim, W_lim] = Newtons_method(rgdx_I{time_n}, vI{time_n}, omegaI{time_n}, kn{time_n},...     %evaluation of residuals and corrections
        knI{time_n}, kapaI{time_n}, gamaI{time_n}, omegaI{time_n+1}, vI{time_n+1}, v_bar{time_n+1}, omega_bar{time_n+1}, ng, n_order,...                  %in the function Newtons_method
        Xg, Wg, ele_c, elem_len, dt, rot_ini, C11, C12, C22, Jp, rho, A, FM_elem, FM_node, n_elem, N_node, fix, dof_free);

        iter = iter + 1;
        fprintf(' -> %2i %17.5e %16.5e \n',iter,W_lim,f_lim);% iteration increment
 
    end                                                                                    % END Loop over Newton iteration 
   
    % Incrementing the quantities after convergence for time t_n+1
    for ele = 1:n_elem                    
        ind_eL = n_order*(ele-1) + 1;
        ind_eR = n_order*ele + 1;
        
        [rg{time_n+1}(:,ind_eL:ind_eR), rgdx{time_n+1}(:,ind_eL:ind_eR), kn{time_n+1}(:,ind_eL:ind_eR), v{time_n+1}(:,ind_eL:ind_eR), omega{time_n+1}(:,ind_eL:ind_eR),...
         qn{time_n+1}(:,ind_eL:ind_eR), gama{time_n+1}(:,ind_eL:ind_eR), kapa{time_n+1}(:,ind_eL:ind_eR)] = time_iteration_update(dt, rot_ini, v{time_n}(:,ind_eL:ind_eR),...
         v_bar{time_n+1}(:,ind_eL:ind_eR), omega{time_n}(:,ind_eL:ind_eR), omega_bar{time_n+1}(:,ind_eL:ind_eR), rg{time_n}(:,ind_eL:ind_eR), rgdx{time_n}(:,ind_eL:ind_eR),...
         kn{time_n}(:,ind_eL:ind_eR), gama{time_n}(:,ind_eL:ind_eR), kapa{time_n}(:,ind_eL:ind_eR), ele_c(ele,:), n_order);

    end
    
    t_i = t_i + dt;                             % time step
    time_n = time_n+1;
    pt(time_n) = t_i;


end
toc

%% Final plots
[f3d, TMe, fxyz] = plots(1, 2, time_n, rg, pt, T_energy, n_node);
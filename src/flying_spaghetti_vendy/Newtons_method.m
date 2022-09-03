function [rgdx_I1, vI1, omegaI1, knI1, kapaI1, gamaI1, v_bar, omega_bar, lim_value, T_energy, f_lim, W_lim] = Newtons_method(rgdx_I0, vI0, omegaI0, kn0, knI0, kapaI0,...
gamaI0, omegaI1, vI1, v_bar, omega_bar, ng, n_order, Xg, Wg, ele_c, elem_len, dt, rot_in, C11, C12, C22, Jp, rho, A, FM_elem, FM_node, n_elem, N_node,...
fix, dof_free)

    %Initializing Global matrices
    K_cons = zeros(6*N_node,6*N_node);                  % Global stiffness matrix
    M_cons = zeros(6*N_node,6*N_node);                  % Global mass matrix
    f_cons = zeros(6*N_node,1);                         % Global residual vector
    K_int = zeros(6*(n_order+1),6*(n_order+1),n_elem);  % Multidimensional array of element stiffness matrix
    M_int = zeros(6*(n_order+1),6*(n_order+1),n_elem);  % Multidimensional array of element mass matrix
    f_int = zeros(6*(n_order+1),n_elem);                % Multidimensional array of element R.H.S vector
    vel_omega = zeros(6*N_node,1);                      % initializing the vector of corrections

    for ele = 1:n_elem                                  % Loop over elements

        ind_eL = n_order*(ele-1) + 1;
        ind_eR = n_order*ele + 1;

        % EVALUATION of TANGENT MATRIX AND VECTOR OF RESIDUALS
        
        [K_int(:,:,ele), M_int(:,:,ele), f_int(:,ele), energy_ele(:,ele)] = Gauss_int(Xg ,Wg, ele_c(ele,:), elem_len(ele,1), elem_len(ele+1,1), dt, vI1(:,:,ele), v_bar(:,ind_eL:ind_eR),...
        vI0(:,:,ele), rgdx_I0(:,:,ele), omegaI1(:,:,ele), omega_bar(:,ind_eL:ind_eR), omegaI0(:,:,ele), gamaI0(:,:,ele), kapaI0(:,:,ele), kn0(:,ind_eL:ind_eR), knI0(:,:,ele),...
        rot_in*C11(2:4,2:4)*rot_in', rot_in*C12(2:4,2:4)*rot_in', rot_in*C22(2:4,2:4)*rot_in', rot_in*Jp*rot_in', A, rho,...
        FM_elem(6*(ele-1)+1:6*(ele+1),1), n_order, ng);
    
        ind = 6*n_order;

        i = 1:6*(n_order+1);

        j = 1:6*(n_order+1);
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
		%  Assembling the global system from data at the element level  %
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        K_cons(ind*(ele-1)+i,ind*(ele-1)+j) = K_cons(ind*(ele-1)+i,ind*(ele-1)+j) + K_int(i,j,ele);
        M_cons(ind*(ele-1)+i,ind*(ele-1)+j) = M_cons(ind*(ele-1)+i,ind*(ele-1)+j) + M_int(i,j,ele);
        f_cons(ind*(ele-1)+i,1) = f_cons(ind*(ele-1)+i,1) + f_int(i,ele);

    end                                                                % Loop over all elements
    T_energy= sum(energy_ele);                                         % Total mechanical energy
    %Rotational transformation of the external moment vector
    for jn = 1:N_node 
        f_ext(6*(jn-1)+1:6*jn,1) = transformed_external_loads(dt, FM_node(6*(jn-1)+1:6*jn,1), omega_bar(:,jn), kn0(:,jn));       
    end
    
    f_cons =  (f_cons- f_ext);                                          % total residual

    % Supported degrees of freedom need to be excluded from the system
    M_cons(:,fix)=[];
    M_cons(fix,:)=[];
    K_cons(:,fix)=[];
    K_cons(fix,:)=[];
    f_cons(fix)=[];
    A_cons = M_cons+K_cons;

    vel_omega(dof_free,1) = -A_cons\f_cons;                           % Vector of corrections
    lim_value = max([norm(vel_omega), norm(f_cons)]);                 % max(norm(residuals),norm(corrections))
    W_lim = norm(vel_omega);
    f_lim = norm(f_cons);
    for k = 1: N_node                                                 % updating velocities (at t_n+1/2) for the next iteration
        vel = vel_omega(6*(k-1)+1:6*k,1);
        v_bar(:,k) = v_bar(:,k) + vel(1:3,1);
        omega_bar(:,k) = omega_bar(:,k) +  vel(4:6,1);      
    end
    
    % updating quantities (positions, velocities, rotations, strains) for the next iteration
    for ele = 1:n_elem
        ind_eL = n_order*(ele-1) + 1;
        ind_eR = n_order*ele + 1;
        [rgdx_I1(:,:,ele), vI1(:,:,ele), omegaI1(:,:,ele), knI1(:,:,ele), kapaI1(:,:,ele), gamaI1(:,:,ele)] = newton_iteration_update(dt, ng, Xg,...
        ele_c(ele,:), elem_len, rgdx_I0(:,:,ele), knI0(:,:,ele), kapaI0(:,:,ele), gamaI0(:,:,ele), omega_bar(:,ind_eL:ind_eR),...
        v_bar(:,ind_eL:ind_eR), vI0(:,:,ele), omegaI0(:,:,ele), ele, n_order);
    end
end
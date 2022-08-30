function [r_n1, rdash_n1, q_n1, v_n1, omega_n1, qn1, gama_n1, kapa_n1] = time_iteration_update(dt, R0, v_n0, v_bar, omega_n0, omega_bar,...
 r_n0, rdash_n0, q_n0, gama_n0, kapa_n0, ele_c, n_order)
    
    
for i = 1: n_order+1
    
    [p,pd] = shape_functions(ele_c,ele_c(i));
    [~, omega_barI, vbdxI, omegabdxI] = Vel_AngVel(n_order, p, pd, v_bar, omega_bar);
    
    [q1, ~] = quaternion_exp(2*dt,omega_barI);
    [q12, q12s] = quaternion_exp(dt,omega_barI);
    
    Q12 = Quaternion_product(q_n0(:,i), q12);
    Q12s = conj_quat(Q12);
    
    t0 = Spurrier(R0);
    q0 = quaternion_rep(t0);
    
    q_n1(:,i) = Quaternion_product(q_n0(:,i),q1);
    qn1(:,i) = Quaternion_product(q0,q_n1(:,i));
    
    r_n1(:,i) = r_n0(:,i) + dt*v_bar(:,i);
    v_n1(:,i) = 2*v_bar(:,i) - v_n0(:,i);
    omega_n1(:,i) = 2*omega_barI - omega_n0(:,i);
    
    rdash_n1(:,i) = rdash_n0(:,i) + dt*vbdxI;
    r_dash_mid = rdash_n0(:,i) + (dt/2)*vbdxI;
    
    T = dq_omega_taylor(dt,omega_barI);
    T_omd = T*[0;omegabdxI];
    
    delgama1 = Quaternion_product(Q12s,Quaternion_product(vbdxI,Q12));
    delgama2 = Quaternion_product(Q12s,Quaternion_product(r_dash_mid,Q12));
    gama_n1(:,i) = gama_n0(:,i) + dt*(delgama1 + skewsymm(delgama2)*[0;omega_barI]);
    
    K12 = Quaternion_product(q12s,Quaternion_product(kapa_n0(:,i),q12)) + 2*Q_left(q12s)*T_omd;
    kapa_n1(:,i) = kapa_n0(:,i) + dt*([0;omegabdxI] - skewsymm(omega_barI)*K12);
    
end


end
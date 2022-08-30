function [rgdxI1, vI1, omegaI1, knI1, kapaI1, gamaI1] = newton_iteration_update(dt, ng, Xg, elec, elem_len, rgdxI0, knI0, kapaI0, gamaI0,...
omega_bar, v_bar, vI0, omegaI0, j, n_order)

h = dt/2;
x1 = elem_len(j,1);
x2 = elem_len(j+1,1);

for i = 1:ng
    
    Xgt(i) = ((x2 - x1)*Xg(i)+(x2 + x1))/2;
    [p, pd] = shape_functions(elec, Xgt(i));    
    [v_barI, omega_barI, vbdxI, omegabdxI] = Vel_AngVel(n_order, p, pd, v_bar, omega_bar);
 
    vI1(:,i) = 2*v_barI - vI0(:,i);
    omegaI1(:,i) = 2*omega_barI - omegaI0(:,i);  
    
    [delq_1,delq_1s] = quaternion_exp(2*dt,omega_barI);
    [delq_12,delq_12s] = quaternion_exp(dt,omega_barI);

    T = dq_omega_taylor(dt,omega_barI);
    T_omd = T*[0;omegabdxI];
    
    knI1(:,i) = Quaternion_product(knI0(:,i),delq_1);
       
    k_12I = Quaternion_product(knI0(:,i),delq_12);
    k_12Is = conj_quat(k_12I);
    
    kapa_12I = Quaternion_product(delq_12s,Quaternion_product(kapaI0(:,i),delq_12)) + 2*Q_left(delq_12s)*T_omd;
    kapaI1(:,i) = kapaI0(:,i) + dt*([0;omegabdxI] - skewsymm(omega_barI)*kapa_12I);
   
    rgdxI1_mid = rgdxI0(:,i) + h*vbdxI;
    rgdxI1(:,i) = rgdxI0(:,i) + dt*vbdxI;
    
    delgama1 = Quaternion_product(k_12Is,Quaternion_product(vbdxI,k_12I));
    delgama2 = Quaternion_product(k_12Is,Quaternion_product(rgdxI1_mid,k_12I));

    gamaI1(:,i) = gamaI0(:,i) + dt*(delgama1 + skewsymm(delgama2)*[0;omega_barI]);
    
end
    

end
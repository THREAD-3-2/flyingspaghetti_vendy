function [V_barI, Omega_barI, V_bar_dash, Omega_bar_dash] = Vel_AngVel(n, p, pd, V_bar, Omega_bar)

V_barI = 0;
Omega_barI = 0;
V_bar_dash = 0;
Omega_bar_dash = 0;

for i = 1:n+1
    
    V_barI = p(i).*V_bar(:,i) + V_barI;
    V_bar_dash = pd(i).*V_bar(:,i) + V_bar_dash;
    
    Omega_barI = p(i).*Omega_bar(:,i) + Omega_barI;
    Omega_bar_dash = pd(i).*Omega_bar(:,i) + Omega_bar_dash;
       
end

end
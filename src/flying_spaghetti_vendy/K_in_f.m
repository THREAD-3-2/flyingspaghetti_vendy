function [Keq, Meq, feq, e_eq]= K_in_f(elec, dt, vI_n1, v_bar, vI_n, rdx_In, omegaI_n1, omega_bar, omegaI_n, gamaI_n, kapaI_n, knI, C_11, C_12, C_22,...
A, rho, Jp, x, j, n)

C11 = eye(4);
C12 = eye(4);
C22 = eye(4);
JP = eye(4);

C11(2:4,2:4) = C_11;
C12(2:4,2:4) = C_12;
C22(2:4,2:4) = C_22;
JP(2:4,2:4) = Jp;
C21 = C12';

%% predictor multiplied with shape functions

[p, pd]  = shape_functions(elec,x);
[v, omega, v_dash, omega_dash] = Vel_AngVel(n, p, pd, v_bar, omega_bar);
[qexp, qexp_s] =  quaternion_exp(dt,omega);

%% Linearization of rotational quaternion

T =  dq_omega_taylor(dt,omega);
T_dash = dxq_omega_taylor(dt,omega,omega_dash);
T_omd = T*[0;omega_dash];

qn12 = Quaternion_product(knI(:,j),qexp);
qn12_s = conj_quat(qn12);

Qmat = Q_right(qn12_s)*Q_left(qn12);
QTmat = Q_left(qn12_s)*Q_right(qn12);
   
   
%% Gama and Kapa calculations (n+1)

rdx_In1 = rdx_In(:,j) + (dt/2)*v_dash;                      
delkapa_n = Quaternion_product(qexp_s,Quaternion_product(kapaI_n(:,j),qexp));
delgama2 =  Quaternion_product(qn12_s,Quaternion_product(rdx_In1,qn12));
delgama1 = Quaternion_product(qn12_s,Quaternion_product(v_dash,qn12));

kapa_barI = delkapa_n + 2*Q_left(qexp_s)*T_omd; 
delq = Q_left(knI(:,j))*T;
    
gama_n1I = gamaI_n(:,j) + dt*(delgama1 + skewsymm(delgama2)*[0;omega]);
kapa_n1I = kapaI_n(:,j) + dt*([0;omega_dash] + skewsymm(kapa_barI)*[0;omega]);

%% n_bar, M_bar

NG = (1/2)*(C11*(gama_n1I + gamaI_n(:,j)) + C12*(kapa_n1I + kapaI_n(:,j)));
MG = (1/2)*(C21*(gama_n1I + gamaI_n(:,j)) + C22*(kapa_n1I + kapaI_n(:,j)));
Ng = Quaternion_product(qn12,Quaternion_product(NG,qn12_s));

%% Skewsymmetric matrices

Sd_gama2 = skewsymm(delgama2);
S_omega = skewsymm(omega);
SKn12 = skewsymm(kapa_barI);
S_MG = skewsymm(MG);
S_NG = skewsymm(NG);
S_JpOm = skewsymm(Jp*omega);

%% Linearized expressions (n + 1/2), Ng, MG

dgama_n12O = (Q_left(delgama2) - Q_right(delgama2))*Q_left(qexp_s)*T;
dgama_n12Vdash = (dt/2)*QTmat;

dgama_n1O = (dt)*((Q_left(delgama1) - Q_right(delgama1))*Q_left(qexp_s)*T + Sd_gama2 - S_omega*dgama_n12O);
dgama_n1Vdash = (dt)*(QTmat - S_omega*dgama_n12Vdash);
    
delKn12_O = (Q_left(delkapa_n) - Q_right(delkapa_n))*Q_left(qexp_s)*T...
 - 2*Q_right(Quaternion_product(qexp_s,T_omd))*Q_left(qexp_s)*T+2*Q_left(qexp_s)*T_dash;
delKn12_Odash = 2*Q_left(qexp_s)*T;
    
delKn1_O = dt*(SKn12 - S_omega*delKn12_O);
delKn1_Odash = dt*(eye(4,4) - S_omega*delKn12_Odash);

NG_Vdash = 1/2*(C11*dgama_n1Vdash);
NG_O = 1/2*(C11*dgama_n1O + C12*delKn1_O);
NG_Odash = 1/2*(C12*delKn1_Odash);
    
Ng_O = (Q_right(Ng) - Q_left(Ng))*Q_right(qn12_s)*delq + Qmat*NG_O;
Ng_Odash = Qmat*NG_Odash;
Ng_Vdash = Qmat*NG_Vdash;
    
gNg_O = Sd_gama2*NG_O - S_NG*dgama_n12O;
gNg_Odash = Sd_gama2*NG_Odash;
gNg_Vdash = Sd_gama2*NG_Vdash - S_NG*dgama_n12Vdash;

MG_O = (1/2)*(C21*dgama_n12O + C22*delKn1_O);
MG_Odash =(1/2)*C22*delKn1_Odash;
MG_Vdash = (1/2)*(C21*dgama_n1Vdash);
  
KMG_O = SKn12*MG_O - S_MG*delKn12_O;
KMG_Odash = SKn12*MG_Odash - S_MG*delKn12_Odash;
KMG_Vdash = SKn12*MG_Vdash;

Jp_O = S_omega*JP - S_JpOm;

%% Energy equation

delWk = (1/2)*(vI_n1(:,j)'*(rho*A*vI_n1(:,j)) + omegaI_n1(:,j)'*Jp*omegaI_n1(:,j));

delWd = (1/2)*(gama_n1I')*C11*gama_n1I + (1/2)*(kapa_n1I')*C22*kapa_n1I;

e_eq = delWk + delWd;
%% Assembly of element matrices
Meq = zeros(6*(n+1),6*(n+1));
Keq = zeros(6*(n+1),6*(n+1));
feq = zeros(6*(n+1),1);

for i = 1:(n+1)
    
    for k = 1:(n+1)
        
        Meq(6*(i-1)+1:6*i-3,6*(k-1)+1:6*k-3) = (2*rho*A*p(i)*p(k)*eye(3,3));
        Meq(6*(i-1)+4:6*i,6*(k-1)+4:6*k) = ((2*Jp)*p(k) + dt*reduce_dim(Jp_O)*p(k))*p(i);
        
        Keq(6*(i-1)+1:6*i-3,6*(k-1)+1:6*k-3) = dt*reduce_dim(Ng_Vdash)*pd(i)*pd(k);
        Keq(6*(i-1)+1:6*i-3,6*(k-1)+4:6*k) = dt*reduce_dim(Ng_O*p(k) + Ng_Odash*pd(k))*pd(i);
        
        Keq(6*(i-1)+4:6*i,6*(k-1)+1:6*k-3) = dt*reduce_dim(MG_Vdash*pd(k) - KMG_Vdash*pd(k)-gNg_Vdash*pd(k))*p(i);
        Keq(6*(i-1)+4:6*i,6*(k-1)+4:6*k) = dt*reduce_dim((MG_O*p(k) + MG_Odash*pd(k))*pd(i) - (KMG_O*p(k) + KMG_Odash*pd(k))*p(i)...
        -(gNg_O*p(k) + gNg_Odash*pd(k))*p(i));
         
    end
    
    feq(6*(i-1)+1:6*i-3,1) = reduce_dim(rho*A*([0;vI_n1(:,j)]-[0;vI_n(:,j)])*p(i) + dt*(Ng*pd(i)));
    feq(6*(i-1)+4:6*i,1) = reduce_dim(JP*([0;omegaI_n1(:,j)]-[0;omegaI_n(:,j)])*p(i)- dt*(SKn12)*MG*p(i) + dt*(S_omega)*(JP*[0;omega])*p(i)+...
    dt*MG*pd(i) - dt*(Sd_gama2*NG)*p(i));
    
end

end
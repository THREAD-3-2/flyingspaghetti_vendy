function [Kint, M_int, f_int, energy_int] = Gauss_int(Xg ,Wg, elec, x1, x2, dt, vI_n1, v_bar, vI_n, rdx_In, omegaI_n1, omega_bar, omegaI_n, gamaI_n, kapaI_n, qn,...
    knI, C11, C12, C22, Jp, A, rho, ext_FM, n_order, ng)

f_ext = [ext_FM(4:6) ext_FM(10:12)];

f_int = 0;
K_int = 0;
M_int = 0;
energy_int = 0;
for j = 1:ng
	Wgtf = Wg(j)*(x2 - x1)/2;
	Xgtf = ((x2 - x1)*Xg(j)+(x2 + x1))/2;

    [Keq, Meq, feq, energy_eq] = K_in_f(elec, dt, vI_n1, v_bar, vI_n, rdx_In, omegaI_n1, omega_bar, omegaI_n, gamaI_n,...
    kapaI_n, knI, C11, C12, C22, A, rho, Jp, Xgtf, j, n_order);
	
    K_int = Wgtf.*Keq + K_int;
	f_int = Wgtf.*feq + f_int;
    M_int = Wgtf.*Meq + M_int;
    energy_int = Wgtf.*energy_eq + energy_int;
end

[m,n] = size(K_int);

qnn(:,1) = qn(:,1);
qnn(:,2) = qn(:,end);

[qexp1,qexp_s1] = quaternion_exp(dt,omega_bar(:,1));
qn121 = Quaternion_product(qnn(:,1), qexp1);
qn12_s1 = conj_quat(qn121);
QT_mat1 = Q_left(qn12_s1)*Q_right(qn121);
T1 = dq_omega_taylor(dt,omega_bar(:,1));

KM_linear{1} = (Q_left(QT_mat1*[0;f_ext(:,1)]) - Q_right(QT_mat1*[0;f_ext(:,1)]))*Q_left(qexp_s1)*T1;

[qexp2,qexp_s2] = quaternion_exp(dt,omega_bar(:,end));
qn122 = Quaternion_product(qnn(:,2), qexp2);
qn12_s2 = conj_quat(qn122);
QT_mat2 = Q_left(qn12_s2)*Q_right(qn122);
T2 = dq_omega_taylor(dt,omega_bar(:,end));

KM_linear{2} = (Q_left(QT_mat2*[0;f_ext(:,2)]) - Q_right(QT_mat2*[0;f_ext(:,2)]))*Q_left(qexp_s2)*T2;

K_ext = zeros(m,n);
K_ext(4:6,4:6) = KM_linear{1}(2:4,2:4);
K_ext(end-2:end,end-2:end) = KM_linear{2}(2:4,2:4);

Kint = K_int - K_ext;

end
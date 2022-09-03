function f_ext = transformed_external_loads(dt, ext_FM, omega, qn)

f_ext = ext_FM;

[qexp, qexp_s] = quaternion_exp(dt,omega);
qn12 = Quaternion_product(qn,qexp);
qn12_s = conj_quat(qn12);
QTmat = Q_right(qn12)*Q_left(qn12_s);
f_ext(4:6,1) = QTmat(2:4,2:4)*ext_FM(4:6,1); 

end
    
    
    
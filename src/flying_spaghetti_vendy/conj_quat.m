function q_conj = conj_quat(q)
q_conj = zeros(4,1);
q_conj(1) = q(1);
q_conj(2:4) = -1*q(2:4);

end
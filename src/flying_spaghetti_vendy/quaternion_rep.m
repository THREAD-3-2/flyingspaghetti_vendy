function Q_theta = quaternion_rep(theta)

theta_norm = (theta'*theta)^0.5;

if theta_norm<eps
    
    theta_norm = eps;

end

factor_1 = cos(theta_norm/2);

factor_2 = (1/theta_norm)*sin(theta_norm/2);

Q_theta = [factor_1;factor_2*theta];

end
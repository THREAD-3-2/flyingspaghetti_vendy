function R = Rodrigues(theta)

% Rotation matrix corresponding to the vector theta

theta_norm = norm(theta);

if theta_norm<eps
    theta_norm = eps;
end

% skew symmetric matrix theta
Theta = reduce_dim(skewsymm(theta));

factor1 = sin(theta_norm)/theta_norm;
factor2 = (1-cos(theta_norm))/(theta_norm^2);

R = eye(3) + factor1*Theta + factor2*Theta*Theta;

end
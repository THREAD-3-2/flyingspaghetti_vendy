function [Q, Qs] = quaternion_exp(dt,omega)

h = dt/2;

onorm = norm(omega);
if onorm > 0.1
    
    fac1 = cos(h*onorm/2);
    fac2 = (omega/onorm)*sin(h*onorm/2);
else   
    fac1 = 1 - ((h^2)*(onorm^2)/8) + ((h^4)*(onorm^4)/384) - (((h*onorm)^6)/46080) + (((h*onorm)^8)/10321920);
    fac2 = ((h/2) - ((h^3)*(onorm^2)/48) + ((h^5)*(onorm^4)/3840) - ((h^7)*(onorm^6)/645120) + (((h^9)*(onorm^8))/185794560)).*omega;
	
end
Q = [fac1;fac2];
Qs = [fac1;-fac2];

end
function T = dq_omega_taylor(dt,omega)

onorm = norm(omega);
T0 = zeros(4,4);
T1 = zeros(4,4);

I43= eye(4);
T0(1,2:4) = (-dt/4)*omega;

T1(2:4,2:4) = omega*omega';
    
if onorm < 0.1
    A0 =((dt/4) - ((dt^3)*(onorm^2)/384) + ((dt^5)*(onorm^4)/122880) - ((dt^7)*(onorm^6)/82575360));
    A1 =( (-(dt^3)/192) + ((dt^5)*(onorm^2)/30720) - ((dt^7)*(onorm^4)/13762560) + ((dt^9)*(onorm^6)/11890851840));        
else    
    A0 = (1/onorm)*sin(dt*onorm/4);
	A1 = (dt/(4*(onorm^2)))*cos(dt*onorm/4) - (1/(onorm^3))*sin(dt*onorm/4);
end    

T = A0*I43 + A0*T0 + A1*T1;

end
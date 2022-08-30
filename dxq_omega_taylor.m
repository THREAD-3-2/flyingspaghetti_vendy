function dq_Tdx = dxq_omega_taylor(dt,omega,omega_dash)

h = dt/2;

I43 = zeros(4);
T0 = zeros(4);
T1 = zeros(4);
T0_dash = zeros(4);
T1_dash = zeros(4);

I43(2:4,2:4)= eye(3);
T0(1,2:4) = omega;
T1(2:4,2:4) = omega*omega';
T0_dash(1,2:4) = omega_dash;
T1_dash(2:4,2:4) = omega*omega_dash' + omega_dash*omega';
onorm = norm(omega);
omomdash = dot(omega,omega_dash);

if onorm > 0.1
	a0 = (1/onorm)*sin(h*onorm/2);
	a1 = (h/(2*(onorm^2)))*cos(h*onorm/2) - (1/(onorm^3))*sin(h*onorm/2);
	b0 = a1*omomdash;
	b1 = -(h/2)*a1*omomdash;
	b2 = -(h/2)*a0;
	b3 = -((h^2)/4)*a0*(omomdash/(onorm^2))-3*a1*(omomdash/(onorm^2));
	b4 = a1;	
else	
	a0 = ((dt/4) - ((dt^3)*(onorm^2)/384) + ((dt^5)*(onorm^4)/122880) - ((dt^7)*(onorm^6)/82575360));
	a1 = ((-(dt^3)/192) + ((dt^5)*(onorm^2)/30720) - ((dt^7)*(onorm^4)/13762560) + ((dt^9)*(onorm^6)/11890851840));
	b0 = a1*omomdash;
	b1 = (-h/2)*a1*omomdash;
	b2 = (-h/2)*a0;
	b3 = (((h^5)/480) - ((h^7)*(onorm^2)/26880) + ((h^9)*(onorm^4)/3870720) - ((h^11)*(onorm^6)/1021870080))*(omomdash);
	b4 = a1;		
end

dq_Tdx = b0*I43 + b1*T0 + b2*T0_dash + b3*T1 + b4*T1_dash;

end


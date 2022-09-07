% Skew symmetric matrix from a vector

function Sm = skewsymm(v)
[m,n] = size(v);

if m==3

	Sm = [0, 0, 0, 0;...
		0, 0, -v(3), v(2);...
		0, v(3), 0, -v(1);...
		0, -v(2), v(1), 0];
else
	
	Sm = [0, 0, 0, 0;...
		0, 0, -v(4), v(3);...
		0, v(4), 0, -v(2);...
		0, -v(3), v(2), 0];
		
end

end
	
   
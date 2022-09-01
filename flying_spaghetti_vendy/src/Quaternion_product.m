function Qp = Quaternion_product(a,b)

if (length(a)== 3) && (length(b) == 4)
   
    r = -(a(1)*b(2) + a(2)*b(3)+a(3)*b(4));
	i = b(1)*a(1) + a(2)*b(4)-b(3)*a(3);
	j = b(1)*a(2)-a(1)*b(4)+b(2)*a(3);
	k = b(1)*a(3)+a(1)*b(3)-b(2)*a(2);
	
	Qp = [r;i;j;k];
    
end

if(length(a)== 4) && (length(b) == 3)
    
    r = -(a(2)*b(1) + a(3)*b(2)+a(4)*b(3));
	i = b(1)*a(1) + a(3)*b(3)-b(2)*a(4);
	j = b(2)*a(1)-a(2)*b(3)+b(1)*a(4);
	k = b(3)*a(1)+a(2)*b(2)-b(1)*a(3);
	
	Qp = [r;i;j;k];
    
end

if (length(a)== 4) && (length(b) == 4)
   
    r = a(1)*b(1)-(a(2)*b(2) + a(3)*b(3)+a(4)*b(4));
	i = a(1)*b(2) + a(2)*b(1)+a(3)*b(4)-b(3)*a(4);
	j = a(1)*b(3)+b(1)*a(3)-b(4)*a(2)+b(2)*a(4);
	k = a(1)*b(4)+b(1)*a(4)+a(2)*b(3)-b(2)*a(3);
	
	Qp = [r;i;j;k];
    
end

if(length(a)== 3) && (length(b) == 3)
    
    r = -(a(1)*b(1) + a(2)*b(2)+a(3)*b(3));
	i = a(2)*b(3)-b(2)*a(3);
	j = b(1)*a(3)-a(1)*b(3);
	k = a(1)*b(2)-b(1)*a(2);
	
	Qp = [r;i;j;k];
    
end

end


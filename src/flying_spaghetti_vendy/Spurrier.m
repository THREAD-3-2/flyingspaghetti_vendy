%% Spurrier Algorithm

% Calculates the rotation vector of theta rotation given by the matrix R
function theta=Spurrier(R)

    tr = trace(R);
    perm=[1 2 3
   	  2 3 1
	  3 1 2];

    [x,i]=max([R(1,1),R(2,2),R(3,3),tr]);

    if (i==4)
        q0=(tr+1)^0.5/2;
        q=1/(4*q0)*[R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)];
    else
        j=perm(i,2);
        k=perm(i,3);
   
        q(i,1)=(R(i,i)/2 + (1-tr)/4)^0.5;
   
        q0=(R(k,j)-R(j,k))/4/q(i);

        q(j,1)=(R(j,i)+R(i,j))/4/q(i);
        q(k,1)=(R(k,i)+R(i,k))/4/q(i);
    end

    norma=norm(q);
    if norma<eps
        norma=1;
    end   

    theta=2*acos(q0)*q/norma;
    
end

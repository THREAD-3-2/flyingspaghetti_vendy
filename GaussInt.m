function [xg,wg]=GaussInt(n)

% evaluates the nodes and weights for Gauss quadrature rule on  [-1,1]
% 
%   n       number of integration points
%   xg      coordinates of integration points (column)
%   wg      weights (row)
%
% integral of a function f on interval [-1,1]  is replaced by the summation sum(wg(i)*f(xg(i))=wg*f(xg) 
% !vectorization is considered here!

if n==1
    wg=2;
    xg=0;
elseif n>=2
	for i=0:n-2
       b(i+1,:)=(i+1)/(((2*i+1)*(2*i+3))^0.5 ); 
	end
	
	K=diag(b,1);
	K=K+diag(b,-1);
	
	[Lve,Vo]=eig(K);
      
	xg=diag(Vo);
	wg=Lve(1,:).^2*2;
	[xg,I]=sort(xg);
	wg=wg(I);  
else
    wg=[];
    xg=[];
    return
end

end
   
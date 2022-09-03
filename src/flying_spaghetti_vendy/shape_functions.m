function [p,pd] = shape_functions(x,X)

n = length(x);
p = zeros(n,1);
pd = zeros(n,1);

for i = 1:n
    
    U = 1;
    L = 1;
    for j = 1:n
        
        if i~=j
            
            U = U*(X-x(j));
            L = L*(x(i)-x(j));
            
        end
        
    end
     p(i) = U/L;
    
end

for j = 1:n
    y = 0;
    for l=1:n
        if not(l==j)
            TEMP = 1/(x(j)-x(l));
            for m=1:n
                if not(m==j) && not(m==l)
                    
                    TEMP = TEMP*(X-x(m))/(x(j)-x(m));
                end
            end
            y = y + TEMP;
        end
    end
    
    pd(j) = y;
    
end

end
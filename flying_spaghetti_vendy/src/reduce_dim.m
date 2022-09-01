function V = reduce_dim(v)


[m,n] = size(v);

if n == m
    
    V = v(2:4,2:4);
    
elseif n<m
    
    V = v(2:4,1);
    
end

end
function h = triangular_load(t)
h = 0;
if t < 2.5  
    h = (200/2.5)*t;   
elseif t < 5   
    h = 400 - ((200/2.5)*t);   
end
end
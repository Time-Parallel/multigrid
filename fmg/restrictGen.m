function R = restrictGen(xf)

nf = length(xf);
nc = nf/2;

R = zeros(nc,nf);

c = 0.5;
o = 0.25;

R(1,1)  = c;
R(1,2)  = o;
R(1,nf) = o;

for i = 2:nc
   
    ind = 2*i-1;
    
    R(i,ind-1) = o;
    R(i,ind) = c;
    R(i,ind+1) = o;
    
end

end
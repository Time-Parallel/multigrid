function P = prolongGen(xc)

nc = length(xc);
nf = 2*nc;

P = zeros(nf,nc);

c = 1;
o = 0.5;
P(1,1) = c;
P(2,1) = o;
for i = 2:nc
   
    ind = 2*i-1;
    
    P(ind-1,i) = o;
    P(ind,i)   = c;
    P(ind+1,i) = o;
    
end
P(nf,1) = o;
end
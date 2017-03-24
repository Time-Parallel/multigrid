clear all;
close all;
clc

N = 288;



L = zeros(N,1);
R = zeros(N,1);
for n = 2:N-1
    
    sl = sum(1:n-1);
    sr = sum(n+1:N);
     L(n) = sl;
     R(n) = sr;
    
    if sl == sr
        n
    end
end

figure
plot(1:N,L,'.r',1:N,R,'.k')


for n = 50:500
   
    t = sqrt(0.5*(n*n+n));
    
    if t == floor(t)
       n
       t
    end
    
    
end
N = 5;

niter = 100000;

rolls = 0;
hits = 0;
bin = rand(N,niter);
kid = rand(N,niter);
for i = 1:niter
    
    rolls = rolls + 1;
    
   
    
    [ibin cbin] = sort(bin(:,i));
    [ikid ckid] = sort(kid(:,i));
    
    if sum(cbin==ckid) == 0
        hits = hits + 1;
    end
end

hits
hits/rolls
   
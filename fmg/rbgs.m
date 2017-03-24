function x = rbgs(lam,b,niter,xstar,w)

N = length(b);

rm = [N 2:2:N-1];
r  = 1:2:N;
rp = 2:2:N;

km = 1:2:N;
k = 2:2:N;
kp = [3:2:N 1];

x = zeros(N,1);

a1 = w;
a2 = 1-w;

for n = 1:niter

    x(r) = a1*x(r) + a2*(b(r) + lam*(x(rm)-x(rp)));
    x(k) = a1*x(k) + a2*(b(k) + lam*(x(km)-x(kp)));
   
    figure(360)
    plot(1:N,xstar,'k',1:N,x,'--r')
    pause(0.1)
    
end
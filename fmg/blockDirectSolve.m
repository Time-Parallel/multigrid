function x = blockDirectSolve(x,b,l,Ax,LXI)

global T
% global a
global CFL
global nu1
% global e
%global dtau_global
global omega
%global CFSO_t
%global DFSO_t
global D_glob
global C_glob
global Nx

n       = 2^(l+1);
dt      = T/n;
dtau    = omega*CFL*dt;

% Direct
%A = genLHS(FSO,a,e,n,dt);
%[xx D] = resGen(DFSO_t,a,e,n,dt,1);
%[Q yy] = resGen(CFSO_t,a,e,n,dt,1);
A = eye(n) + C_glob{l} + D_glob{l};
L = eye(n) + dtau*A;
LI = L\eye(n);
dx = zeros(size(x));
dxx = dx;
res = dx;
%res =  b-A*x;
%l20 =  norm(res,2)/sqrt(n);
for i = 1:nu1
    for j = 2:Nx
        res(j,1:n) =  b(j,1:n)-(A*x(j,1:n)')';
    end
        %l2 = norm(res,2)/sqrt(n);
    
    for j = 2:Nx
        dx(j,1:n) = (LI*(dtau*res(j,1:n)'))';
    end
    
    for k = 1:n
        dxx(1:Nx,k) = LXI*dx(1:Nx,k);
    end
    
    x = x + dxx;
    %disp(['Residual on level ' num2str(l) ' : ' num2str(l2/l20) ])
end



end
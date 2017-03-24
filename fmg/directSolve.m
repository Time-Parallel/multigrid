function x = directSolve(x,b,l)

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
%res =  b-A*x;
%l20 =  norm(res,2)/sqrt(n);
for i = 1:nu1
    res =  b-A*x;
    %l2 = norm(res,2)/sqrt(n);
    dx = LI*(dtau*res);
    x = x + dx;
    %disp(['Residual on level ' num2str(l) ' : ' num2str(l2/l20) ])
end
% Jacobi Iteration

% im      = [n 1:n-1];
% ij      = 1:n;
% ip      = [2:n 1];
% i2dt    = 1/2/dt;
% 
% for i = 1:nu1
%     x(ip)
%     x(im)
%     b(ij)
%     
%     v = dtau*( b(ij)-i2dt*(x(ip)-x(im)) );
%     
%     x = v;
% end

end


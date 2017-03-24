function x = presmooth(x,b,l)

global T
global a
global CFL
global nu1
global e
global FSO

n       = 2^(l+1);
dt      = T/n;
dtau    = CFL*dt/a;

% Direct
%A = genLHS(FSO,a,e,n,dt);
[xx D] = resGen(4,a,e,n,dt);
[Q yy] = resGen(2,a,e,n,dt);
A = Q + D;
L = eye(n) + dtau*A;
LI = L\eye(n);
disp(['New Smoother '])
for i = 1:nu1
    res =  b-A*x;
    l2 = norm(res,2)/sqrt(n);
    dx = LI*(dtau*res);
    x = x + dx;
    disp(['Residual on level ' num2str(l) ' : ' num2str(l2) ])
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


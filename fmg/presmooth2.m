function x = presmooth2(x,b,l)

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
% L = eye(n) + dtau*A;
% LI = L\eye(n);
disp(['New Smoother '])
alpha = [1/4 1/6   3/8 1/2     1]; 
beta  = [  1   0 14/25   0 11/25];
for i = 1:nu1
    res =  b-A*x;
    l2 = norm(res,2)/sqrt(n);
    %dx = LI*(dtau*res);
    %x = x + dx;
    un = x;
    q0 = Q*un-b;
    d0 = D*un;
    
    u1 = un - alpha(1)*dtau*(q0+d0);
    
    q1 = Q*u1-b;
    d1 = D*u1;
    
    u2 = un - alpha(2)*dtau*(q1 + beta(2)*d1 + (1-beta(2))*d0);
    
    q2 = Q*u2-b;
    d2 = D*u2;
    
    u3 = un - alpha(3)*dtau*(q2 + beta(3)*d2 + (1-beta(3))*d1);
    
    q3 = Q*u3-b;
    d3 = D*u3;
    
    u4 = un - alpha(4)*dtau*(q3 + beta(4)*d3 + (1-beta(4))*d2);
    
    q4 = Q*u4-b;
    d4 = D*u4;
    
    u5 = un - alpha(5)*dtau*(q4 + beta(5)*d4 + (1-beta(5))*d3);
    
    x = u5;
    
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


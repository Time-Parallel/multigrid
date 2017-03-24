function x = postsmooth(x,b,l)

global T
global a
global CFL
global nu2
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
disp(['New Post Smoother '])
for i = 1:nu2
    res =  b-A*x;
    l2 = norm(res,2)/sqrt(n);
    dx = LI*(dtau*res);
    x = x + dx;
    disp(['Residual on level ' num2str(l) ' : ' num2str(l2) ])
end


end


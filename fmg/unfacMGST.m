function x = unfacMGST(x,b,l,Ax,LXI)

global T
%global FSO
%global a
global e
global lmin
%global ux
global nu1
global nu2
global dtau_global
global direct_smooth
global direct_solve
% global CFSO_t
global DFSO_t
% global C_glob
% global D_glob
% global P_glob
% global R_glob
global Nx


% figure(31)
% plot(1:length(x),x,'r',1:length(x),ux,'--k')
 
if l > lmin
 
    x = unfacSmoothMGST(x,b,l,nu2,Ax,LXI);
  
    
    % Restrict Residual
    nl  = 2^(l+1);
    dt  = T/nl;
    n   = nl;
%     A   = eye(nl) + C_glob{l} + D_glob{l};

    idt = 1/dt;

    ca  =  0.5*idt*dtau_global*[0   -1   0   1  0];

    if DFSO_t == 2
        cs  = e*idt*dtau_global*[0  -1   2  -1  0];
    elseif DFSO_t == 4
        cs  = e*idt*dtau_global*[1  -4   6  -4  1];
    else
       disp('Incompatible Dissipation Operator'); 
    end
    ca(3) = 1 + ca(3);
    
    c = ca + cs;  
    
    jmm = [n-1 n 1:n-2];
    jm  = [n 1:n-1];
    j   = 1:n;
    jp  = [2:n 1];
    jpp = [3:n 1 2];

    res(:,j) = b(:,j)-(c(1)*x(:,jmm) + c(2)*x(:,jm) + c(3)*x(:,j) + c(4)*x(:,jp) + c(5)*x(:,jpp));
    
    for k = 1:n
        res(1:Nx,k) = res(1:Nx,k) - dtau_global*Ax*x(1:Nx,k);
    end
    
    %disp(['Level: ' num2str(l) ' :: L2: ' num2str(norm(res,2)/sqrt(Nx*n))])
    
    r   = blockRestrict(res);
    [dummy, nl1] = size(r);
    
    v   = zeros(Nx,nl1);
    v   = unfacMGST(v,r,l-1,Ax,LXI);
    
    % Prolongation
    p   = blockProlong(v);
    x   = x + p;
    
    % Postsmoothing
    x = unfacSmoothMGST(x,b,l,nu2,Ax,LXI);
else
    x = unfacSmoothMGST(x,b,l,nu1,Ax,LXI);
end

end

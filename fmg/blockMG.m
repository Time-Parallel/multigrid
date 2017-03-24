function x = blockMG(x,b,l)

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
    if direct_smooth == 1
        for ii = 2:Nx 
            duu = x(ii,:)';
            duu = directSolve(duu,b(ii,:)',l);
            x(ii,:) = duu;
        end
    else
        x = blockSmoothMG(x,b,l,nu2);
    end
    
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
   
  %  norm(res,2)/sqrt(Nx*n)
    
    r   = blockRestrict(res);
    [dummy, nl1] = size(r);
    
    v   = zeros(Nx,nl1);
    v   = blockMG(v,r,l-1);
    
    % Prolongation
    p   = blockProlong(v);
    x   = x + p;
    
    % Postsmoothing
    if direct_smooth == 1
         for ii = 2:Nx 
            duu = x(ii,:)';
            duu = directSolve(duu,b(ii,:)',l);
            x(ii,:) = duu;
        end
    else
        x = blockSmoothMG(x,b,l,nu2);
    end
else
    if direct_solve == 1
        for ii = 2:Nx 
            duu = x(ii,:)';
            duu = directSolve(duu,b(ii,:)',l);
            x(ii,:) = duu;
        end
        
    else
        x = blockSmoothMG(x,b,l,nu1);
    end
end

end

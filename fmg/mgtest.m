function x = mgtest(x,b,l)

%global T
%global FSO
%global a
%global e
global lmin
%global ux
global nu1
global nu2
%global dtau_global
global direct_smooth
global direct_solve
%global CFSO_t
%global DFSO_t
global C_glob
global D_glob
global P_glob
global R_glob

% figure(31)
% plot(1:length(x),x,'r',1:length(x),ux,'--k')
 
if l > lmin
    if direct_smooth == 1
        x = directSolve(x,b,l);
    else
        x = smoothMG(x,b,l,nu2);
    end
    % Restrict Residual
    nl  = 2^(l+1);
    % dt = T/nl;

    A   = eye(nl) + C_glob{l} + D_glob{l};
    res = b-A*x;
    r   = R_glob{l}*res;
    nl1 = length(r);
    
    v   = zeros(nl1,1);
    v   = mgtest(v,r,l-1);
    
    % Prolongation
    p   = P_glob{l-1}*v;
    x   = x + p;
    
    % Postsmoothing
    if direct_smooth == 1
        x = directSolve(x,b,l);
    else
        x = smoothMG(x,b,l,nu2);
    end
else
    if direct_solve == 1
        x = directSolve(x,b,l);
    else
        x = smoothMG(x,b,l,nu1);
    end
end

end

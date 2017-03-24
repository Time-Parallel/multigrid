function [U,L2,cyc] = iterate(U,l,cyc,ncyc,tol,L2)

%% Global Variables
global w
global T    % Time Period
global a
global e_t
global e
global CFSO_t
global DFSO_t
global FSO_x
global CFL 
global Nx
global dx
global nplot
global tsolve
global mgst
global nmgcyc
global x
global dtau_global
global l2r0
global Nl

N = 2^(l+1);

% Temporal Inputs
t       = linspace(0,T-T/N,N)';
dt      = T/N;
dt0     = dt;

% Spatial Inputs
e_x     = e;

%% Initialization
XX      = zeros(Nx,N);
TT      = XX;
UX      = XX;
DUDX    = UX;
B       = UX;
RES     = UX;
DU      = UX;
DUU     = UX;
for j = 1:Nx
    for k = 1:N
        XX(j,k)     = x(j);
        TT(j,k)     = t(k);
        UX(j,k)     = sin(w*(x(j)-a*t(k)));
        DUDX(j,k)   = -w*a*cos(w*(x(j)-a*t(k)));
        B(j,k)      = -DUDX(j,k);
    end
end

U(1,:) = UX(1,:);

figure(l)
pcolor(XX,TT,U)
shading flat
title(['U0 at Level ' num2str(l) ])

% Temporal Operators
[xx, D] = resGen(DFSO_t,1,e_t,N,dt,1);  % Dissipative
[C, yy] = resGen(CFSO_t,1,e_t,N,dt,1);  % Convective
A = C + D;


dtau = CFL/(1/dt + a/dx);
dtau_global = dtau;

% Spatial Operators
[Cx, Dx] = resGen(FSO_x,a,e_x,Nx,dx,0);
Ax = Cx + Dx;

err = abs(U-UX);
l2e = norm(err,2)/sqrt(N*Nx);

for j = 2:Nx
    RES(j,1:N) = RES(j,1:N) + (A*U(j,1:N)')';
end
for k = 1:N
   RES(:,k) = RES(:,k) + Ax*U(:,k); 
end

l2r = norm(RES,2)/sqrt(N*Nx);
%l2r0=l2r;



LX  = eye(Nx) + dtau*Ax;
LXI = LX\eye(Nx);

LT  = eye(N) + dtau*A;
LTI = LT\eye(N);

tic
dtau0 = dtau;

while l2r > tol && cyc < ncyc
    cyc         = cyc + 1;

    % Evaluate RHS
    for j = 2:Nx
        RES(j,1:N)  = (A*U(j,1:N)')';
    end
    for k = 1:N
        RES(:,k)    = RES(:,k) + Ax*U(:,k); 
    end    
    
    if mgst == 0
        % Invert Temporal Part
        if tsolve == 1
            for j = 2:Nx
                rhs = -dtau*RES(j,1:N)';
                duu = LTI*rhs;
                DUU(j,1:N) = duu(1:N);
            end
        elseif tsolve == 2
           for j = 2:Nx
              duu   = 0*U(j,1:N)';
              b     = -dtau*RES(j,1:N)';

              for mgcyc = 1:nmgcyc
                duu = mgtest(duu,b,l);
              end

              DUU(j,1:N) = duu(1:N);
           end
        elseif tsolve == 3
            for j = 2:Nx
                rhs = -dtau*RES(j,1:N)';
                duu = LT\rhs;
                DUU(j,1:N) = duu(1:N);
            end
        elseif tsolve == 4
            duu = 0*U;
            b = -dtau*RES;
            for mgcyc = 1:nmgcyc
                duu = blockMG(duu,b,l);
            end
            DUU(2:Nx,1:N) = duu(2:Nx,1:N);
        elseif tsolve == 5
            for j = 2:Nx
                rhs = -dtau*RES(j,1:N)';
                %duu = thomasPeriodic(LT,rhs,5,j);
                duu = thomasPeriodic(LT,rhs,API2,uu,vv,2,j);
                DUU(j,1:N) = duu(1:N);
            end
        end

        % Invert Spatial Part
        for k = 1:N
            rhs = DUU(1:Nx,k);
            du = LXI*rhs;
            DU(1:Nx,k) = du(1:Nx);
        end
    elseif mgst == 1
       b = -dtau*RES;
       for mgcyc = 1:nmgcyc
        DU = blockMGST(DU,b,l,Ax,LXI); 
      %  disp(norm(DU,2)/sqrt(N*Nx))
       end
    elseif mgst == 2
        b = -dtau*RES;
        DUU = 0*DU;
        for mgcyc = 1:nmgcyc
           DU = unfacMGST(DU,b,l,Ax,LXI);
        %   disp(norm(abs(DU-DUU),2)/sqrt(Nx*N))
        %   DUU = DU;
        end
    end
    
    % Update
    U(2:Nx,:) = U(2:Nx,:) + DU(2:Nx,:);
    
    % Plottng
    if mod(cyc,nplot)==0
        figure(33)
        subplot(6,2,[1 3])
        pcolor(XX,TT,U)
        shading flat
        colorbar
        title(cyc)
        
        subplot(6,2,[5 7])
        pcolor(XX,TT,log10(abs(U-UX)))
        shading flat
        colorbar
        title('log_{10}(e)')
        
        subplot(6,2,[9 11])
        pcolor(XX,TT,log10(abs(RES)))
        shading flat
        colorbar
        title('log_{10}(R)')
        
        subplot(6,2,[2 4 6])
        semilogy(1:cyc,L2(1:cyc,1))
        xlim([1 cyc+3*nplot])
        ylabel('|e|_2')
        
        subplot(6,2,[8 10 12])
        semilogy(1:cyc,L2(1:cyc,2))
        xlim([1 cyc+3*nplot])
        ylim([10^(-16) 10^2])
        ylabel('|R|_2')
    end
    
    err         = abs(U-UX);
    l2e         = norm(err,2)/sqrt(N*Nx);
    l2r         = norm(RES,2)/sqrt(N*Nx)/l2r0;
    disp(['Iteration: ' num2str(cyc) ' :: Residual: ' num2str(l2r) ' :: Error: ' num2str(l2e)])
    L2(cyc,1) = l2e;
    L2(cyc,2) = l2r;
end

clear all;
close all;
clc

%% Global Variables
global N
global w
global T    % Time Period
global nu1  % Presmooth
global nu2  % Postmooth
global a
global e_t
global e
global CFSO_t
global DFSO_t
global FSO_x
global CFL 
global lmin
global direct_smooth
global direct_solve
global dtau_global
global omega
global smooth_factor
global C_glob
global D_glob
global P_glob
global R_glob
global Nx
global dx
global dt0
global nplot
global x
global tsolve
global mgst
global nmgcyc
global mglevel
global l2r0
global alpha
global beta
global Nl

%% Inputs & Parameters

% Global Inputs
tol     = 1e-12;    % Convergence Tolerance
a       = 1;        % Convective Speed
e       = 1/32;    % Artificial Dissipation Coefficient (need for time and space)
nu1     = 1;        % Number of smoothing iterations on coarsest grid
nu2     = 1;        % Number of smoothing iterations on finer grids
ncyc    = 10000;    % Maximum number of mg cycles
CFL     = 3;       % Space-Time CFL
nplot   = 1000;     % Plot solution every nplot cycles

% Temporal Inputs

% Solver Type:
% 1) Direct Solve (No MG)
% 2) 
% 3) 
% 4) 
% 5) 
% 6) 
% 7) 
tsolve  = 4;    
N       = 1024;     % Temporal Degrees of Freedom (must be power of two now but will modify for any number)
w       = 2*pi;     % Frequency
T       = 2*pi/w;   % Temporal Period
e_t     = e;        % Temporal artificial dissipation coefficient
CFSO_t  = 2;        % Convective Order of Accuracy (only 2nd now - hardcoded as 2nd in lower routines)
DFSO_t  = 4;        % AD Order of Accuracy (2nd and fourth now - hardcoded as above)
mglevel = 3;        % How many levels of MG
mgst    = 1;         
nmgcyc  = 1;        % How many MG cycles per (nonlinear/global) iteration
fmgcyc = 200;       % Number of MG cycles on coarse grids to initialize
fmgtol = 1e-6;      % Convergence tolerance on coarser grids
coeffs = 1;         % Standard/Optimized (see Kaveh Husseini) Multidstage Coefficients

if coeffs == 1
    alpha = [1/4 1/6   3/8 1/2     1]; 
    beta  = [  1   0 14/25   0 11/25];
elseif coeffs == 2
    alpha = [0.2312 0.1679 0.3678 0.4996 1]; 
    beta  = [1      0      0.5413 0      0.3222];
end

omega   = 100;      % Timestep multiplier for direct solver
t       = linspace(0,T-T/N,N)'; % time
Nl      = floor(log2(N)) - (1 + (CFSO_t-2)/2); % Overall number of MG levels
l       = Nl;       
nl      = 2^(l+1);
dt      = T/nl;     % Physical timestep
dt0     = dt;       
lmin    = max(1,Nl-mglevel+1);% Go down to 3 points for 2nd order, 5 points for fourth order, etc
lstart  = lmin; % Start on level lstart for FMG
direct_solve    = 0; % direct solves on coarsest (not implemented for all schemes yet)
direct_smooth   = 0; % direct solve as smoother on finer grids (not implemented for all schemes yet)
smooth_factor   = 3.5/CFL; % Like a CFL for the multistage scheme (still confirming if scaled correctly)

% Spatial Inputs
Nx      = 101;      % Spatial Degrees of Freedom
Lx      = 0.5;     % Spatial Domain size
e_x     = e;       % Spatial AD
FSO_x   = 6;       % Spatial Order of Accuracy (2,4,6)
x       = linspace(0,Lx,Nx)';
dx      = x(2)-x(1);

%% Initialization
XX      = zeros(Nx,N);
TT      = XX;
UX      = XX;
U       = UX;
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
        U(j,k)      = 0;
        DUDX(j,k)   = -w*a*cos(w*(x(j)-a*t(k)));
        B(j,k)      = -DUDX(j,k);
    end
end

U(1,:) = UX(1,:);

figure
pcolor(XX,TT,UX)
shading flat
title('U0')

% Temporal Operators
[xx, D] = resGen(DFSO_t,1,e_t,N,dt,1);  % Dissipative
[C, yy] = resGen(CFSO_t,1,e_t,N,dt,1);  % Convective
A = C + D;

dtau = CFL/(1/dt + a/dx);
%dtau = CFL/(a/dx);
dtau_global = dtau;

for level = lmin:Nl
    
    nn = 2^(level+1);
    ndt = T/nn;
    nt  = linspace(0,T-T/nn,nn)';
    
    % Dissipation Operator
    [xx, D] = resGen(DFSO_t,1,e_t,nn,ndt,1);  % Dissipative
    D_glob{level} = dtau*D;
    
    % Convective Operator
    [C, yy] = resGen(CFSO_t,1,e_t,nn,ndt,1);  % Convective
    C_glob{level} = dtau*C;
    
    % Prolongation Operator
    P = prolongGen(nt);
    P_glob{level} = P;
    
    % Restriction Operator
    R = restrictGen(nt);
    R_glob{level} = R;

end

% Spatial Operators
[Cx, Dx] = resGen(FSO_x,a,e_x,Nx,dx,0);
Ax = Cx + Dx;

cyc = 0;
err = abs(U-UX);
l2e = norm(err,2)/sqrt(N*Nx);

for j = 2:Nx
    RES(j,1:N) = RES(j,1:N) + (A*U(j,1:N)')';
end
for k = 1:N
   RES(:,k) = RES(:,k) + Ax*U(:,k); 
end

l2r = norm(RES,2)/sqrt(N*Nx);
l2r0 = l2r;
l2r=l2r/l2r0;

L2 = zeros(10000,2);

LX  = eye(Nx) + dtau*Ax;
LXI = LX\eye(Nx);

LT  = eye(N) + dtau*A;
LTI = LT\eye(N);


tic
dtau0 = dtau;

%% Get Solution to Coarsest Level
uf = U;

if lstart ~= Nl
    for l = Nl:-1:lstart+1 

        uc = blockRestrict(uf);
        clear uf;
        uf = uc;

    end
    U = uc;
end

lcyc = zeros(size(lstart:Nl));

for l = lstart:Nl
        
    if l == Nl;
        mgtol = tol;
        fmgcyc = ncyc;
    else
        mgtol = fmgtol;
    end
    
    scyc = cyc;
    
    %lmin = max(lmin,l-2)
    
    [U,L2,cyc] = iterate(U,l,cyc,cyc+fmgcyc,mgtol,L2);

    lcyc(l) = cyc-scyc+1;
    
    if l ~= Nl
        uf = blockProlong(U);
        clear U;
        U = uf;
    end
    
end

t_exec = toc
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

figure
surf(XX,TT,U)
shading flat

figure
subplot(2,1,1)
semilogy(L2(:,1))

subplot(2,1,2)
semilogy(L2(:,2))

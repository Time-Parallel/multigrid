clear all;
close all;
clc

N = 256;
e = 1/8;
DFSO = 4;
CFSO = 2;

w = 2*pi;
T = 2*pi/w;
tt = linspace(0,T-T/N,N)';

dt = T/(N-1);

L = 1;
Nx = 50;

dx = L/(Nx-1);

x = rand(1);

CFL = 2;

dtau = CFL*dt;

[xx, D] = resGen(DFSO,1,e,N,dt,1);  % Dissipative
[C, yy] = resGen(CFSO,1,e,N,dt,1);  % Convective

A = C+D;

u = zeros(N,1);
a = 1;
xb = rand(1);
b = -dtau*(w*cos(w*(xb-a*tt))+(C+D)*u);
ux = sin(w*tt);
dudt = w*cos(w*tt);

maxiter = 100;

% for n = 1:100
% 

% 
%     b = -dtau*(dudx + (C+D)*u);
%A = rand(r  = zeros(N,maxiter+1);

tol = 1e-14;
nplot = 1;
c=0;
for n = 1:1

    r = zeros(N,maxiter+1);
    rh = r;
    v  = r;
    p  = r;
    %b = rand(N,1);

    % A = [1 2 3;4 7 9; 8 7 2];
    % b = [1;2;3];

    %xstar = A\b;
    
    x = zeros(N,1);
    
    b = dudt;
    
    r(:,1)  = b - A*x;
    rh(:,1) = r(:,1);

    rho   = 1;
    alpha = 1;
    omega = 1;

    v(:,1) = 0;
    p(:,1) = 0;
    
    for i = 1:maxiter
        c = c + 1;
        rhop = rho;
        rho = rh(:,1)'*r(:,i);

        beta = (rho/rhop)*alpha/omega;

        p(:,i+1) = r(:,i) + beta*(p(:,i)-omega*v(:,i));

        v(:,i+1) = A*p(:,i+1);

        alpha = rho/(rh(:,1)'*v(:,i+1));

        h = x + alpha*p(:,i+1);



        s = r(:,i)-alpha*v(:,i+1);

        t = A*s;

        omega = (t'*s)/(t'*t);
        xp = x;
        x = h + omega*s;

        res = b-A*x;
        l2r = norm(res,2)/sqrt(N);
        L2(c) = l2r;

        if l2r < tol
            i
            break
        end
        r(:,i+1) = s-omega*t;

    end
    u = x;
    %u = u + x;
    if mod(n,nplot)==0
        figure(33)
        plot(tt,ux,'ok',tt,u,'.r')
        ylim([-1.2 1.2])
    end
end
%     
figure
semilogy(abs(L2))
%     u = u + x;
%     
%     figure(33)
%     plot(tt,ux,'ok',tt,x,'.-r')
%     title(n)

% end
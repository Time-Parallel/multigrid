clear all;
close all;
clc

N = 64;

A = zeros(N);
u = zeros(N,1);
v = u;


w = 2*pi;
T = 2*pi/w;
t = linspace(0,T-T/N,N)';
dt = T/N;

CFL = 0.5;

dtau = CFL*dt;


c = 0;
b = 1 + dtau/dt;
a = -dtau/dt;

b = 1;
a = -dtau/2/dt;
c = dtau/2/dt;


A(1,1) = b;
A(1,N) = a;

for i = 2:N-1
    A(i,i)   = b;
    A(i,i-1) = a;
    A(i,i+1) = c;
end

A(N,N) = b;
A(N,1) = c;

u(1) = -b;
u(N) = c;

v(1) = 1;
v(N) = -a/b;

AP = A-u*v';


AI = A\eye(N);
API = AP\eye(N);
API2 = zeros(N);

tol = 0.01;

for i = 1:N
    for j = 1:N
    
        if abs(API(i,j))>tol
            API2(i,j) = API(i,j);
        end
        
    end
end

I1 = API2*AP;
I2 = AP*API2;

% figure
% subplot(3,1,1)
% spy(API)
% subplot(3,1,2)
% spy(API2)
% subplot(3,1,3)
% spy(I1)

figure
spy(API2)

b = rand(N,1);

xstar = AP\b;

x = zeros(N,1);


ncomm =  sum(abs(API2(:,floor(N/2)))>tol);

tol_conv = 1e-14;

res0 = norm(x-xstar,2);
set_conv = 0;
niter = 100;
n_conv = niter;
res = zeros(niter,1);

for n = 1:niter
    
    x = x - API2*(AP*x-b);
    
    res(n) = norm(x-xstar,2)/res0;
    
    if res(n) < tol_conv && set_conv == 0
        set_conv = 1;
        n_conv = n;
    end
     
end

if set_conv == 0
   n_conv = Inf;
end

ops = n_conv*ncomm;

figure
semilogy(res)
title(ops)

G = eye(N) - API2*AP;
    
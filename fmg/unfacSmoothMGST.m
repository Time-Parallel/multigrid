function x = unfacSmoothMGST(x,b,l,nu,Ax,LXI)

global N
global Nx
global T
%global a
global CFL
global e
global dtau_global
global smooth_factor
% global CFSO_t
global DFSO_t
global dx
global a
% global C_glob
% global D_glob

% if nu < 1
%     return
% else

%disp(['Level ' num2str(l) ' :: sweeps = ' num2str(nu) ]);

n       = 2^(l+1);
dt      = T/n;
%dtau_old    = smooth_factor*CFL*dt;
dtau    = smooth_factor*CFL/(dtau_global/dt);

%dtau = smooth_factor*CFL/(1+1/dt+a/dx);
%dtau = CFL/(1/dt+a/dx);
%dtau = 3.5/(1+dtau_global/dt);
%dtau    = smooth_factor*CFL/(dtau_global/dt);
%  D = D_glob{l};
%  C = eye(n) + C_glob{l};

alpha = [1/4 1/6   3/8 1/2     1]; 
beta  = [  1   0 14/25   0 11/25];

idt = 1/dt;

ca  =  0.5*idt*dtau_global*[   -1   0   1   ];

if DFSO_t == 2
    cs  = e*idt*dtau_global*[0  -1   2  -1  0];
elseif DFSO_t == 4
    cs  = e*idt*dtau_global*[1  -4   6  -4  1];
else
   disp('Incompatible Dissipation Operator'); 
end
%ca(2) = 1 + ca(2);

jmm = [n-1 n 1:n-2];
jm  = [n 1:n-1];
j   = 1:n;
jp  = [2:n 1];
jpp = [3:n 1 2];
q   = zeros(Nx,n);
d0  = q;
d2  = q;
d4  = q;
u1 = q;
u2 = q;
u3 = q;
u4 = q;
%disp(['Size of x is ' num2str(size(x)) ' and size of q is ' num2str(size(q))])

% for k = 1:n
%     b(1:Nx,k) = LXI*b(1:Nx,k);
% end
Ax = dtau_global*Ax;
LXI_1 = (eye(Nx) + dtau*alpha(1)*(Ax+eye(Nx)))\eye(Nx);
LXI_2 = (eye(Nx) + dtau*alpha(2)*(Ax+eye(Nx)))\eye(Nx);
LXI_3 = (eye(Nx) + dtau*alpha(3)*(Ax+eye(Nx)))\eye(Nx);
LXI_4 = (eye(Nx) + dtau*alpha(4)*(Ax+eye(Nx)))\eye(Nx);
LXI_5 = (eye(Nx) + dtau*alpha(5)*(Ax+eye(Nx)))\eye(Nx);
for i = 1:nu

    un = x;

    q(:,j)  =                   ca(1)*un(:,jm) + ca(2)*un(:,j) + ca(3)*un(:,jp) - b(:,j);
    d0(:,j) = cs(1)*un(:,jmm) + cs(2)*un(:,jm) + cs(3)*un(:,j) + cs(4)*un(:,jp)+cs(5)*un(:,jpp);
    for k = 1:n
        u1(1:Nx,k)  =   LXI_1*(un(1:Nx,k) - alpha(1)*dtau*(q(1:Nx,k) + d0(1:Nx,k)));
    end
    

    q(:,j)  =                   ca(1)*u1(:,jm) + ca(2)*u1(:,j) + ca(3)*u1(:,jp) - b(:,j);
    for k = 1:n
        u2(1:Nx,k)  =   LXI_2*(un(1:Nx,k) - alpha(2)*dtau*(q(1:Nx,k) + d0(1:Nx,k)));
    end

    q(:,j)  =                   ca(1)*u2(:,jm) + ca(2)*u2(:,j) + ca(3)*u2(:,jp) - b(:,j);
    d2(:,j) = cs(1)*u2(:,jmm) + cs(2)*u2(:,jm) + cs(3)*u2(:,j) + cs(4)*u2(:,jp)+cs(5)*u2(:,jpp);
    rd20    = beta(3)*d2 + (1-beta(3))*d0;
    for k = 1:n
        u3(1:Nx,k)  =   LXI_3*(un(1:Nx,k) - alpha(3)*dtau*(q(1:Nx,k) + rd20(1:Nx,k)));
    end 

    q(:,j)  =                  ca(1)*u3(:,jm) + ca(2)*u3(:,j) + ca(3)*u3(:,jp) - b(:,j);
    for k = 1:n
       u4(1:Nx,k)   =   LXI_4*(un(1:Nx,k) - alpha(4)*dtau*(q(1:Nx,k) + rd20(1:Nx,k)));
    end

    q(:,j)  =                   ca(1)*u4(:,jm) + ca(2)*u4(:,j) + ca(3)*u4(:,jp) - b(:,j);
    d4(:,j) = cs(1)*u4(:,jmm) + cs(2)*u4(:,jm) + cs(3)*u4(:,j) + cs(4)*u4(:,jp) + cs(5)*u4(:,jpp);
    
    for k = 1:n
        x(1:Nx,k)   =   LXI_5*(un(1:Nx,k) - alpha(5)*dtau*(q(1:Nx,k) + beta(5)*d4(1:Nx,k) + (1-beta(5))*rd20(1:Nx,k)));
    end
end

end
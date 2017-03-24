clear all
close all
clc

Nx = 101;
N  = 1024;
Niter = 100;

XF = rand(Nx,N);
XF0 = XF;
xf = XF(1,:)';

ts = zeros(Niter,1);
tb = ts;
tt = ts;
for n = 1:Niter
    tic
    for i = 1:Nx     
    %xc = restrict(XF0(i,1:N)');
        xc = restrict(xf);
        xf = prolong(xc);
    end
    titer = toc;
    ts(n) = titer;
end
%tsweep = toc;
for n = 1:Niter
    tic
    XC = blockRestrict(XF);
    XF = blockProlong(XC);
    titer = toc;
    tb(n) = titer;
end

clear XF
XF = XF0';
for n = 1:Niter
    tic
    XC = blockRestrict(XF);
    XF = blockProlong(XC);
    titer = toc;
    tt(n) = titer;
end

figure
plot(1:Niter,ts,1:Niter,tb,1:Niter,tt)
legend('Sweep','Row','Column')
title(sum(ts)/sum(tb))

figure
plot(1:Niter,tb,1:Niter,tt)
legend('Row','Column')
title(sum(tt)/sum(tb))

%tblock = toc;

%tsweep/Niter
%tblock/Niter
%tsweep/tblock


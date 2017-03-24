function xc = blockRestrict(xf)

[nx, nf] = size(xf);

nc = nf/2;

jc = 0:nc-1;

jm = mod(2*jc-1,nf);
j  = 2*jc;
jp = mod(2*jc+1,nf);

xc(1:nx,jc+1) = 0.25*xf(1:nx,jm+1) + 0.5*xf(1:nx,j+1) + 0.25*xf(1:nx,jp+1);

end
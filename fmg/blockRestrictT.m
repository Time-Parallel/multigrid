function xc = blockRestrictT(xf)

[nf,nx] = size(xf);

nc = nf/2;

jc = 0:nc-1;

jm = mod(2*jc-1,nf);
j  = 2*jc;
jp = mod(2*jc+1,nf);

inx = 1:nx;

xc(jc+1,inx) = 0.25*xf(jm+1,inx) + 0.5*xf(j+1,inx) + 0.25*xf(jp+1,inx);

end
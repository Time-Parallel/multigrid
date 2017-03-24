function xf = blockProlong(xc)

[nx,nc] = size(xc);

nf = 2*nc;
xf(1:nx,1:2:nf) = xc(1:nx,1:nc);
xf(1:nx,2:2:nf) = 0.5*( xc(1:nx,1:nc) + xc(1:nx,[2:nc 1]) ); 

end
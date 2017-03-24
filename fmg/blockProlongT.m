function xf = blockProlongT(xc)

[nc,nx] = size(xc);

nf = 2*nc;

% Dirct Injection
inx = 1:nx;
xf(1:2:nf,inx) = xc(1:nc,inx);
xf(2:2:nf,inx) = 0.5*( xc(1:nc,inx) + xc([2:nc 1],inx) ); 

end
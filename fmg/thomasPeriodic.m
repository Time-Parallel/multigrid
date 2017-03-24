function x = thomasPeriodic(A,b,API2,u,v,niter,iter)

N = length(b);

% y = zeros(N,1);
% q = API2*u;
% x = y;
% for n = 1:niter
% 
%     x = x - (y - ((v'*y)/(1+v'*q))*q);    
% 
%     if n < niter
%         r = A*x-b;
%         y = API2*r;
%     end
% 
% end

x = zeros(N,1);
UVT = u*v';
for n = 1:niter
    x = API2*(b-UVT*x);
    
%     if iter == 10
%     disp(norm(b - A*x,2)/sqrt(N));
%     end
    
end
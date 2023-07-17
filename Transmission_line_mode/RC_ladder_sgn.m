function [f,B,C] = RC_ladder_sgn(N)
%
% System: xdot = f(x) + B*u;
%         y = C x
%
% Input: * N: order of original model

A =-2*spdiags([ones(N,1);0],0,N,N)+spdiags(ones(N,1),1,N,N)+spdiags([ones(N-1,1);0],-1,N,N);
n = @(x) sign(x).*x.^2;

f = @(x) A*x-n(x);

B = zeros(N,1);
B(1) = 1;

C = zeros(1,N);
C(1) = 1;

end
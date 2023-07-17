function fJac = RC_ladder_sgn_Jac(x,N)
%analytical calculation of jacobi matrix of RC_Ladder

g = @(x) -2*x.^2.*dirac(x)-2*x.*sign(x)-ones(N,1)*2;
fJac = spdiags(g(x),0,N,N) + spdiags(ones(N,N),-1,N,N)+ spdiags(ones(N,N),1,N,N);
%fJac = diag(g(x))+ full(spdiags(ones(N,N),-1,N,N)+ spdiags(ones(N,N),1,N,N));
end
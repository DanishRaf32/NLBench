function g = NHT_nonlinearities_spd(x,a,N,A,l)
q = zeros(N, 1);
for i=1:1:length(a)-1 %a1=a(2)
    q=q+a(i+1)*x.^(i+1)/(i+1);
end

G=sparse(N,N);
% G=-spdiags(ones(N,N),-1,N,N)+2*spdiags(ones(N,N),0,N,N)-spdiags(ones(N,N),1,N,N);
G=-spdiags(ones(N,1),-1,N,N)+2*spdiags(ones(N,1),0,N,N)-spdiags(ones(N,1),1,N,N);
G(1,1)=1;

g=G*q;
g=-1*g*A/l;
end
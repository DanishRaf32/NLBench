function fJac = NHT_Jac_spd(x,N,K,E,A,a,l)
% calculation of Jacobian "without" for loops
qJac=sparse(N,N);

for i=1:1:length(a)-1 %a1=a(2)
    qJac=qJac+a(i+1)*x.^i;
end

G=sparse(N,N);
G=-spdiags(ones(N,1),-1,N,N)+2*spdiags(ones(N,1),0,N,N)-spdiags(ones(N,1),1,N,N);
G(1,1)=1;

diagqJac=spdiags(speye(N,N)*qJac*speye(N,1),0,N,N); 

gJac=(A/l)*G*diagqJac;

fJac=-1*gJac;% jacobian of nonlinearity
end

%% OHNE x im Aufruf (similar to implementation 3)?!?
% function fJac = NHT_Jac_spd(N,A,a,l) 
% 
% fJac = @(x) (A/l)*gJac; % jacobian of nonlinearity
%
% end


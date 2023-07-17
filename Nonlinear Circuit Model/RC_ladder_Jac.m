function fJac = RC_ladder_Jac(x,N)
%analytical calculation of jacobi matrix of RC_Ladder

A0 = sparse(N,N);
A0(1,1)=1; 

A1 = spdiags(ones(N-1,1),-1,N,N)-speye(N);
A1(1,1)=0;

A2 = spdiags([ones(N-1,1);0],0,N,N)-spdiags(ones(N,1),1,N,N);

e = @(x) exp(40.0*x);

P = spdiags(e(A0*x),0,N,N);
Q = spdiags(e(A1*x),0,N,N);
R = spdiags(e(A2*x),0,N,N);

J_p = -P+transpose(transpose(A1)*Q-transpose(A2)*R);

fJac = -A0+A1-A2 + 40*J_p;

end


% if loop
%     e1 = e(A1*x);
%     e2 = e(A2*x);
%     Jcurr=sparse(N,N);
%     for i=1:1:N
%         for j=1:1:N
%             Jcurr(j,i)=A1(j,i)*e1(j,1)-A2(j,i)*e2(j,1);
%         end
%     end
%     fJac = A1-A2 + 40*Jcurr;
% 
% else
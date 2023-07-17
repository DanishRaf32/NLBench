function [f,B,C] = RC_ladder(N)
%
% System: xdot = f(x) + B*u;
%         y = C x
%
% Input: * N: order of original model

g = @(x) exp(40.0*x)+x-1;

% A0
A0 = sparse(N,N);
A0(1,1)=1; 
% A1
A1 = spdiags(ones(N-1,1),-1,N,N)-speye(N);
A1(1,1)=0;
% A2
A2 = spdiags([ones(N-1,1);0],0,N,N)-spdiags(ones(N,1),1,N,N);

f = @(x) -g(A0*x)+g(A1*x)-g(A2*x);

B = zeros(N,1);
B(1) = 1;

C = zeros(1,N);
C(1) = 1;

end

%---------------------- Adapted MOR-Wiki Code ----------------------------%
% function [xdot, y] = RC_ladder(N)
% %
% % System: xdot@(x,u)
% %         y@(x)
% %
% %         xdot = @(x,u) -g(A0*x)+g(A1*x)-g(A2*x) + [u;sparse(N-1,1)];
% %         y = @(x) x(1);
% %
% % Input: * N: order of original model
% 
% g = @(x) exp(40.0*x)+x-1;
% 
% A0 = sparse(N,N); A0(1,1)=1; 
% A1 = spdiags(ones(N-1,1),-1,N,N)-speye(N); A1(1,1)=0;
% A2 = spdiags([ones(N-1,1);0],0,N,N)-spdiags(ones(N,1),1,N,N);
% 
% xdot = @(x,u) -g(A0*x)+g(A1*x)-g(A2*x) + [u;sparse(N-1,1)];
% y = @(x) x(1);
% 
% end

%%
%%------------------------ AUXILIARY FUNCTIONS --------------------------%%

% function fJac = RC_ladder_Jac(x,A1,A2)
% %analytical calculation of jacobi matrix
% 
% e = @(x) exp(40.0*x);
% 
% e1 = diag(e(A1*x));
% e2 = diag(e(A2*x));
% 
% J_p = transpose(A1)*e1-transpose(A2)*e2;
%  
% fJac = A1-A2+40*J_p;
% 
% end

% if loop
%     e  = @(x) exp(40.0*x);
%     e1 = e(A1*x);
%     e2 = e(A2*x);
% 
%     Jcurr=sparse(N,N);
%     for i=1:1:N
%         for j=1:1:N
%             Jcurr(j,i)=A1(j,i)*e1(j,1)-A2(j,i)*e2(j,1);
%         end
%     end
%     fJac = A1-A2 + 40*Jcurr;
% end

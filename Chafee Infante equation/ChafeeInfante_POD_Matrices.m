function [E,A,B,C] = ChafeeInfante_POD_Matrices(k,L)

% Chafee Infante benchmark
% E*xdot = A*x + B*u(t) - fcubic(x)

h = L/(k+1);
E = speye(k);
A = sparse(k,k);
B = sparse(k,1);
C = sparse(1,k);

% left boundary
A(1,1) = 1-2/h^2;
A(1,2) = 1/h^2;

% right boundary 
A(k,k-1) = 1/h^2;
A(k,k) = 1-1/h^2;

for i = 2:k-1
   A(i,i-1) = 1/h^2;
   A(i,i) = 1-2/h^2;
   A(i,i+1) = 1/h^2;
end

B(1,1) = 1/h^2;
C(1,k) = 1;
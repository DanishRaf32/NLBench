function [E,A,B,C] = Fitz_Matrices_POD(k,l)
% FitzHugh-Nagumo system 
%
h = l/(k-1);
% E = speye(2*k);
% A = sparse(2*k,2*k);
% B = sparse(2*k,2);

E = speye(3*k);
A = sparse(3*k,3*k);
B = sparse(3*k,2);

%Alternative?
%e = 1;
%R = 0.01;
%gamma = 1;
%g = 0*0.05; 


e = 0.015;
R = 0.5;
gamma = 2;
g = 0.05;

% for i = 1:(k/2)
%    E(i,i) = e;
%    E(i+2*k,i+2*k) = e;
% end

for i = 1:k
   E(i,i) = e;
   E(i+2*k,i+2*k) = e;
end

% left boundary
A(1,1) = -e^2/h^2-0.1;
A(1,2) = e^2/h^2;
A(1,k+1) = -1;
%A(1,2*k+1) = 1.1;%---------------------------------------------------

A(k+1,1) = R;
A(k+1,k+1) = -gamma;

A(2*k+1,1) = 2*g;
A(2*k+1,2*k+1) = -2*e^2/h^2-0.2;

% right boundary
A(k,k-1) = e^2/h^2;
A(k,k) = -e^2/h^2-0.1;
A(k,2*k) = -1;
%A(k,3*k) = 1.1;-----------------------------------------

A(2*k,k) = R;
A(2*k,2*k) = -gamma;

A(3*k,3*k) = -2*e^2/h^2-0.2;
A(3*k,k) = 2*g;

% inner points
for i = 2:k-1
   A(i,i-1) = e^2/h^2;
   A(i,i) = -2*e^2/h^2-0.1;
   A(i,i+1) = e^2/h^2;
   A(i,i+k) = -1;
  % A(i,i+2*k) = 1.1;---------------------------------
   
   A(i+k,i) = R;
   A(i+k,i+k) = -gamma;
   
   A(i+2*k,i) = 2*g;
   A(i+2*k,i+2*k) = -4*e^2/h^2-0.2;
end

B(1,1) = e^2/h;
% B(:,2) = g;
B(1:2*k,2) = g;

C = sparse(2,3*k);
C(1,1) = 1;
C(2,1+k) = 1;
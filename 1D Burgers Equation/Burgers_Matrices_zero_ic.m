function [E,A,H,N,B,C] = Burgers_Matrices_zero_ic(k,mu)
% Viscous Burgers equation
% u_t + u * du/dx = mu * d^2u/dx^2
% subject to zero initial condition u_0(x) = 0
% and boundary control on the left boundary: u(0)=u(t)
% homog. Neumann boundary condition on the right: u'(1)=0
% Semidiscretization with k points leads to:
% du_1/dt = -u_1 * (u_1-u_0)/h + mu * (u_2 - 2*u_1 +u_0)/h^2
% du_i/dt = -u_i * (u_i-u_(i-1))/h + mu * (u_(i+1) - 2*u_i + u_(i-1))/h^2
% du_k/dt = -u_k * (u_k-u_(k-1))/h + mu * (u_(k+1) - 2*u_k + u_(k-1))/h^2
% Measured output is assumed to be the right boundary, i.e. u_k
h = 1/(k+1);
E = speye(k);
A = sparse(k,k);
H = sparse(k,k^2);
N = sparse(k,k);
B = sparse(k,1);
C = sparse(1,k);

% left boundary
A(1,1) = -2*mu/h^2;
A(1,2) = mu/h^2;
% right boundary
A(k,k-1) = mu/h^2;
A(k,k) = -mu/h^2;
% inner points
for i = 2:k-1
    A(i,i-1) = mu/h^2;
    A(i,i) = -2*mu/h^2;
    A(i,i+1) = mu/h^2;
end

% left boundary
H(1,1) = -1/h;
% right boundary
H(k,k^2) = -1/h;
H(k,(k-2)*k+k) = 1/(2*h);
H(k,(k-1)*k+k-1) = 1/(2*h);
% inner points
for i = 2:k
    H(i,(i-1)*k+i) = -1/h;
    H(i,(i-1)*k+i-1) = 1/(2*h);
    H(i,(i-2)*k+i) = 1/(2*h);
end

N(1,1) = 1/h;
B(1) = mu/h^2; 
C(end) = 1;
% C=1/k*ones(1,k);
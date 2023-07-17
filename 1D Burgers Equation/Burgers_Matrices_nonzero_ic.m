function [E,A,H,N,B,C,x0] = Burgers_Matrices_nonzero_ic(k,mu)

% Viscous Burgers equation
% u_t + u * du/dx = mu * d^2u/dx^2
% subject to non-zero initial condition u_0(x) = u_0 and
% homog. Neumann boundary condition on the right: u'(1)=0
% Semidiscretization with k points leads to:
% du_1/dt = mu * (u_2 - u_1)/h^2
% du_i/dt = -u_i * (u_i-u_(i-1))/h + mu * (u_(i+1) - 2*u_i + u_(i-1))/h^2
% du_k/dt = -u_k * (u_k-u_(k-1))/h + mu * (u_(k+1) - 2*u_k + u_(k-1))/h^2
% Measured output is assumed to be the average value of u
h = 1/(k+1);
E = speye(k);
A = sparse(k,k);
H = sparse(k,k^2);
N = sparse(k,k);

I = speye(k);
 x0 = zeros(k,1);
%x0=sign(linspace(0,5,k));
 for i = 1:k
    % x0(i)=1+(sin((2*i/k+1)*pi));
    x0(i)=1-(sin((2*pi/k)*i));
    
 end
% left boundary
A(1,1) = -mu/h^2;
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

C = 1/k*ones(1,k);
B = A*x0 + H*kron(x0,x0);
A = A + H*(kron(I,x0)+kron(x0,I));
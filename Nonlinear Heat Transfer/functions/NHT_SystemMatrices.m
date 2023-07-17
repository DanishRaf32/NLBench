function [E,K,B,C] = NHT_SystemMatrices(N,a,l,A,cp,rho)
%N elements --> N-nodes: right end of the beam is fixed at ambient
%                        temperature 0
%calculatin of system matrices: E x' + K x = B u + g(T)
%with x=T
%a = [a0 a1 ... an] the coefficients of the polynomial representing the
%                   nonlinear heat conductivity kappa

%% calculation of M, B and K Matrix
% K-Matrix is the linear part of the nonlinear stiffness matrix
% ("Nonlinear Heat Transfer Modelling and Reduction", equation 9b)

% =========================================================================
% assembly of matrices
% =========================================================================    
B=sparse(N,2);
M=sparse(N,N);
K=sparse(N,N);


M=spdiags(ones(N,1),-1,N,N)*1/6+2/3*speye(N,N)+spdiags(ones(N,1),1,N,N)*1/6;
M(1,1)=1/3; %overwrite wrong value at position 1,1 with correct one: 1/3

K=spdiags(ones(N,1),-1,N,N)*-1+2*speye(N,N)+spdiags(ones(N,1),1,N,N)*-1;
K(1,1)=1;
K=a(1)*A/l*K;


B(1,1)=A;
B(1,2)=(A*l)/2;
B(2:N,2)=A*l;
% =========================================================================
% Boundary condition: T=0 at x=L
% =========================================================================

% =========================================================================

E=rho*cp*A*l*M;

C=sparse(2,N); % see "NLHeatTransfer p.4"
C(1,1)=1;
C(2,round(N/2))=1;
end


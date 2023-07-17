function [J] = MemsJac(xold,consts,Operators)
% function [J] = MemsJac(xold,consts,Operators)
%
% consts = [S0 E In0 rho0 lambda p0 z0 width dy viscosity dx height b0 del M N];
% Operators = [{INT} {X} {D2} {grady} {X} {X} {gradx2} {gradx} {LAP} {D4}]
%      |f1|
%  f = |f2|
%      |f3|
%
%     |J11 J12 J13|
%  J =|J21 J22 J23|
%     |J31 J32 J33|

M = consts(15); 
N = consts(16);
%order = (M+2)*N; 
S0 = consts(1); 
E = consts(2); 
In0 = consts(3); 
rho0 = consts(4); 
lambda = consts(5); 
p0 = consts(6); 
z0 = consts(7); 
width=consts(8); 
%dy = consts(9); 
viscosity = consts(10);
height = consts(12); 
%dx = consts(11); 
del = consts(14);
rho = rho0*height*width;
%S = S0*height*width;
In = In0*width*height^3;

INT = Operators{1}; 
D4 = Operators{10}; 
D2 = Operators{3}; 
grady = Operators{4}; 
gradx = Operators{8}; 
LAP = Operators{9}; 
gradx2 = Operators{7};

x1 = xold(1:N); x2 = xold(N+1:2*N); x3 = xold(2*N+1:N*(M+2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J11 = -2/3*diag(x2./(x1.^3))*1/del; J11=sparse(J11);
J12 = 1/3*diag(1./(x1.^2))*1/del;   J12=sparse(J12);
J13 = zeros(N,N*M);                 J13=sparse(J13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J21 = -2*diag(x2.^2./(x1.^4))*1/del + 6/rho*diag(x1)*diag( INT*(x3-p0) + S0*height*width*(D2*(x1-z0)) - E*In*(D4*(x1-z0)) )*del...
    + 3/rho*diag(x1.^2)*(S0*height*width*D2 - E*In*D4)*del;
J22 = 4/3*diag((x2/del)./(x1.^3));
J23 = 3/rho*diag(x1.^2)*INT*del;
J21=sparse(J21); J22=sparse(J22); J23=sparse(J23);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J31a = diag(x3)*WeirdDiag( diag((x2/del)./(x1.^4)),M,N) ;
J31b = diag(x3.*(gradx*(x3-p0)))*WeirdDiag((diag(gradx2*x1-z0)+diag(x1)*gradx2+4*lambda*gradx2)/(4*viscosity),M,N);
J31c = diag(x3.*(LAP*(x3-p0)+(gradx*(x3-p0)).^2+(grady*(x3-p0)).^2))*WeirdDiag( diag((2*x1+6*lambda)/(12*viscosity)),M,N);
J31 = J31a + J31b + J31c;
J31=sparse(J31);
%%%
J32 = -diag(x3)*WeirdDiag( diag(1./(3*x1.^3)),M,N)*1/del;
J32=sparse(J32);
%%%
arg331b = (x1+4*lambda)./(4*viscosity).*(gradx2*(x1-z0));
J331 = -DiagHat(((x2/del)./(3*x1.^3)),M,N) + DiagHat(arg331b,M,N)*(diag(gradx*(x3-p0))) + diag(x3)*gradx ;
%%%
J332 = DiagHat(((x1.^2+6*lambda*x1)/(12*viscosity)),M,N)*(diag(LAP*(x3-p0)) + diag(x3)*LAP + 2*diag(gradx*(x3-p0))*gradx+2*diag(grady*(x3-p0))*grady);
%%%
J33 = J331 + J332;
J33=sparse(J33);
%%%%%%%%%%%%%%%
J = [J11 J12 J13 ; J21 J22 J23; J31 J32 J33];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Jtemp = DiagHat(x,M,N)
Jtemp = diag(reshape(repmat(x',M,1),N*M,1));

function Jtemp2 = WeirdDiag(A,M,N)
Jtemp2 = reshape(repmat(A(:)',M,1),M*N,N);

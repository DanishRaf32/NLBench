function [f B] = MemsFunc(xold,consts,Operators)
% function [f B] = MemsFunc(xold,consts,Operators)
%
% consts = [S0 E In0 rho0 lambda p0 z0 width dy viscosity dx height b0 del M N];
% Operators = [{INT} {X} {D2} {grady} {X} {X} {gradx2} {gradx} {LAP} {D4}]
%
%      |f1|
%  f = |f2|
%      |f3|

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
B = Operators{5};

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

x1 = xold(1:N); 
x2 = xold(N+1:2*N); 
x3 = xold(2*N+1:N*(M+2));

f1 = x2./(3*x1.^2)*1/del;

f2 = (2*x2.^2)./(3*x1.^3)*1/del + (3*x1.^2)/(rho).*( INT*(x3-p0)+S0*height*width*(D2*(x1-z0))-E*In*(D4*(x1-z0)))*del;

arg31 = (x2./(3*x1.^3))*1/del;
f31 = -DiagHat((arg31),M,N) * x3;

arg32 = (x1+4*lambda)./(4*viscosity).*(gradx2*(x1-z0));
f32 = DiagHat((arg32),M,N) * (x3.*(gradx*(x3-p0)));

arg33 = ((x1.^2+6*lambda*x1)/(12*viscosity));
f33 = DiagHat((arg33),M,N) * (x3.*(LAP*(x3-p0))+(gradx*(x3-p0)).^2+(grady*(x3-p0)).^2);

f3 = f31 + f32 + f33;

f = [f1 ; f2 ; f3];

%%%%%%%%%%%%%%%%%%%%%%%%%
function Jtemp = DiagHat(x,M,N)
Jtemp = diag(reshape(repmat(x',M,1),N*M,1));

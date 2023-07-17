%% initialize and prepare
clear
clc
close all

%% parameters
N=13;                    %Number of elements
h=10e-6;
w=15e-6;
e0=8.854e-12;           %epsilon_0
% s=200e-9;               %distance of upper and lower beam in undeformed configuration
s=200e-6;
L=100e-6;               %length of beam
rho=2330;
A=150e-12;
I=1.25e-21;
E=1.31e11;
d1=0;
d2=1e-6;

%calculated parameters
n=N+1;
l=L/N;                  %lenght of beam element

%--------------------------------------------------------------------------
%Boundary condition
%--------------------------------------------------------------------------
% Bndry=[1 n];         %-> two sides: y1=yn=0   -> [1 n]
Bndry=1;               %-> fixed left side: y1=theta1=0; ->1

%--------------------------------------------------------------------------
% system order
%--------------------------------------------------------------------------
%   -> order can be 0, 1 or 2: 
%   (0 for calculation of dae with daetoolbox)
order = 1;

%% system
[E1st, A1st, b1st, c1st] = ESB_SystemMatrices(n,N,rho,E,I,l,h,w,A,s,e0,d1,d2,Bndry,order);

P = sparse(4*n,4*n);
g = sparse(4*n,1);

f1st = @(x) ESB_nonlinearities(n,s,Bndry,order,x,e0,l);

%%
% u = @(t) 100*(t>=0.1);
u = @(t) 100*sin(2*pi*5*t);

left = @(t,x) [-A1st E1st]*x - b1st*u(t) - f1st;

%%
PJac=sparse(4*n, 4*n);
gJac=sparse(4*n,4*n);
J = @(t, x) ESB_Jac(n,s,Bndry,order,x,e0,l,PJac,gJac);

%% Simulation
dt = 0.001;
tEnd = 0.3;
tSim = 0:dt:tEnd;

%%
x0 = zeros(size(E1st,1),1);

%% simulation

% Uniform step size h
h = tSim(2)-tSim(1);

% Initial time step and y
tnew = tSim(1) + h;
xcurr = x0;

% Preallocation of y
x = zeros(length(tSim),length(x0));
x(1,:) = x0;

%Numerical or analytical Jacobian:
FJacTimeDep = @(t, x) E1st/h - A1st - J(t,x);
    
tic
% Implicit Euler method
for i = 1:length(tSim)-1
    sys = @(xnew) E1st*((xnew - xcurr)/h) - A1st*xnew - b1st*u(tnew) - f1st(xnew);
    FJac = @(x) FJacTimeDep(tnew, x);
    
    xnew = NewtonRaphson(sys,xcurr,FJac);
    x(i+1,:) = xnew;
    tnew = tnew + h;
    xcurr = xnew;
end
simTime = toc
y=c1st*x.';

figure;
plot(tSim, y);
xlabel('Time (s)');
ylabel('Output y');

X=0:1/N:1; %x-coordinate/L
if Bndry==[1 n]
xSpatial=zeros(n,1);
    for i=1:1:n-length(Bndry)
        xSpatial(i+1) = x(length(x)-3,2*i);    
    end
end
if Bndry==1
    for i=1:1:(n-2)
        if i==1
             xSpatial(i+2) = x(length(x)-3,i);    
        end
        xSpatial(i+2) = x(length(x)-3,2*i-1);    
    end
end

figure;
plot(X, xSpatial);

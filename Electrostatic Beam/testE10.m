%% Benchmark collection of large-scale nonlinear state-space models for dimensionality reduction
% Authors: Danish Rafiq, M. A. Bazaz 
clear
clc
close all

%% Prepare Model
%Load Model
sys_name = 'E10';

%% Load nonlinear model
[M, E, K, F, f2nd, B, C, G, x0, v0] = ElectrostatBeam10;
n2nd = size(M,1);

%% Convert model from second order to first order
E1st = [speye(n2nd,n2nd) sparse(n2nd,n2nd); sparse(n2nd,n2nd) M];
A1st = [sparse(n2nd,n2nd) speye(n2nd,n2nd); -K -E];
f1st = @(x) [sparse(n2nd,1); f2nd(x)];
B1st = [sparse(n2nd,1); B];
C1st = [C, sparse(1,n2nd)];

n1st = size(A1st,1);


%--------------------------------------------------------------------------
% parameters
%--------------------------------------------------------------------------
N=9;                   %Number of elements
h=10e-6;
w=15e-6;
e0=8.854e-12;           %epsilon_0
s=200e-9;               %distance of upper and lower beam in undeformed configuration
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

%Boundary condition 
Bndry=[1 n];            %-> two sides: y1=yn=0   -> [1 n]

% system order
order = 1;              % can be 1 or 2; for simulation order must be 1

% Jacobian of nonlinearity
J = @(t, x) ESB_Jac(n,s,Bndry,order,x,e0,l);

%% Simulation 
% Simulation parameters
dt = 0.01; %Step width
tEnd = 1;
tSim = 0:dt:tEnd;

%Excitation signal
% u = @(t) 1*(t>=0);
u = @(t) 100*sin(2*pi*5*t);

% Initial state vector
x02nd = x0;
x01st = [x0; v0];

fTest2 = f2nd(x02nd);
fTest1 = f1st(x01st);

% Simulation of the nonlinear system
xdot = @(t,x) A1st*x + f1st(x) + B1st*u(t);   
Jac_xdot = @(t,x) A1st + J(t,x);
options = odeset('Mass',E1st,'Jac',Jac_xdot);

% [t_nonlin,x_nonlin] = ode15s(xdot,tSim,x01st,options); %ode15s only for Index-1 DAEs
% [t_nonlin,x_nonlin] = ode23t(xdot,tSim,x01st,options); %ode23t only for Index-1 DAEs

tic
[t,x] = ode23s(xdot,tSim,x01st,options);
simTime = toc
y=C1st*x.';

%--------------------------------------------------------------------------
% plotting the results
%--------------------------------------------------------------------------
% Output y
figure;
plot(tSim, y);
xlabel('Time [s]');
ylabel('Output y [m]');

% displacement y(x,t=tEnd)
X=0:1/N:1; %x-coordinate/L
xSpatial=zeros(n,1);
for i=1:1:n-length(Bndry)
    xSpatial(i+1) = x(length(x)-3,2*i);    
end

figure;
plot(X, xSpatial);
xlabel('x/L [-]');
ylabel('displacement y(x) [m]');

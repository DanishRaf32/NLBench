%% initialize and prepare
clear
clc
close all

%% system

%--------------------------------------------------------------------------
% parameters
%--------------------------------------------------------------------------
N=10;                   %Number of elements
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
%Bndry=1;                %-> fixed left side: y1=theta1=0; ->1
% system order
order = 1;              % can be 1 or 2; for simulation order must be 1

%--------------------------------------------------------------------------
% system
%--------------------------------------------------------------------------
% matrices
Opts.transf2nd = 'I'; %-K
[E1st, A1st, b1st, c1st] = ESB_SystemMatrices(n,N,rho,E,I,l,h,w,A,s,e0,d1,d2,Bndry,order,Opts);

% nonlinearity
f1st = @(x) ESB_nonlinearities(n,s,Bndry,order,x,e0,l);

% Jacobian of nonlinearity
J = @(t, x) ESB_Jac(n,s,Bndry,order,x,e0,l);

%% Simulation
%--------------------------------------------------------------------------
% simulation parameters
%--------------------------------------------------------------------------
dt = 0.001;
tEnd = 0.3;
tSim = 0:dt:tEnd;

x0 = zeros(size(E1st,1),1);

% excitation signal
u = @(t) 100*sin(2*pi*5*t);
%u=@(t) 100;
% Simulation of the nonlinear system
xdot = @(t,x) A1st*x + f1st(x) + b1st*u(t);   
Jac_xdot = @(t,x) A1st + J(t,x);
options = odeset('Mass',E1st,'RelTol',1e-6,'Jac',Jac_xdot);

% [t_nonlin,x_nonlin] = ode15s(xdot,tSim,x0,options); %ode15s only for Index-1 DAEs
% [t_nonlin,x_nonlin] = ode23t(xdot,tSim,x0,options); %ode23t only for Index-1 DAEs

tic
[t,x] = ode23s(xdot,tSim,x0,options);
simTime = toc
y=c1st*x.';

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
xlabel('x/L [-]');
ylabel('displacement y(x) [m]');


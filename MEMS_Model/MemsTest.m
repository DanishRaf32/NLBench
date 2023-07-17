%% Benchmark collection of large-scale nonlinear state-space models for dimensionality reduction
% Authors: Danish Rafiq, M. A. Bazaz 
clear all; clc; close all;
%Spatial discretization
N = 30;  % # of points in Length discretization
M = 15;  % # of points in Width discretization
% Initial conditions
x(:,1)=rand((2+M)*N,1);
% random vector to linearize around
%x0test=rand((2+M)*N,1);
x0(1:N,1)= 2.3;
x0(N+1:2*N,1)=0;
x0(2*N+1:(2+M)*N,1)=0.1013;
x0test=x0;
% define constants and matrix finite-difference operators
[consts, operators] = MemsConstsAndOperators(N,M);

% compute function evaluation 
[f B] = MemsFunc(x0test,consts,operators);

% compute Jacobian
[J] = MemsJac(x0test,consts,operators);
 
%%
omega=2*pi*10e3;
%omega=(pi/3)*4e4;
tf=5e-4;%final time
dt=1e-6;
tspan=0:dt:tf-dt;
k=round(tf/dt);
%u=@(t) (14*cos(omega*t)).^2;
%u=@(t) (8.58*sin(60*pi*t)).^(2);
%u=49*ones(k,1);
%u=9.^2;
%u=@(t) (7*cos(4*pi*t)).^2;
u= @(t) (7*cos((4*pi/30)*4e4*t)).^2;
%u=(7*cos((4*pi/30)*4e4*time)).^2;
%u= @(t)(3*cos(2*omega*t)+7*cos(0.5*omega*t)).^2;
%u=@(t) ((7*cos(2*omega*t))+(3*cos(0.5*omega*t))).^2;

xdot=@(t,x) MemsFunc(x,consts,operators)+B*u(t);
options.solver='NewtonRaphson';
tic
[tsol,xsol]=ode15s(xdot,tspan,x0test,options);
Sim_time=toc
xsol=xsol';
C=zeros(1,size(B,1));
C(1,N/2)=1;
y=C*xsol;
plot(tsol,y,'k','Linewidth',3)
xlabel('time ($s$)','Interpreter','LaTeX');
ylabel('Beam deflection at center-point (microns)','Interpreter','LaTeX');
set(gca,'FontSize',15,'TickLabelInterpreter','LaTeX')
grid on


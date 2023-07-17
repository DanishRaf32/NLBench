% Author: Danish Rafiq
% This MATLAB script simulates 1D Burger's model for different viscosity
% values
% 1D Burgers' PDE: Ut + U*Ux = nu * Uxx  for x =[0 L], t=[0 tEnd]
% A single control input u(t) is applied on left boundary
% Parameters Range : [nu_range1,nu_range2]
clear all; close all; clc;
%% simulation parameter
k = 500;             % No. of spatial points
L = 5;               % Length of spatial variable
x=linspace(0,L,k);   % Spatial variable
t0 = 0;              % Initial time
tEnd = 15;           % Final time 
dt = 0.01;           % Time Discretization Step
tSim = t0:dt:tEnd;   % Time span vector
nu_1=0.01;           % Initial Parameter Range vlaue
nu_2=5;             % Final Parameter Range value
nu_range=6;          % Number of Parameter intervals
i=1;                 % Global index
u = @(t) exp(-t);    % Test input 
for nu=linspace(nu_1,nu_2,nu_range)
disp(['simulating BE for Viscosity nu=' num2str(nu)])
[E,A,H,N,B,C] = Burgers_Matrices_zero_ic(k,nu); % E=I
nQB = size(A,1);
xdotQB = @(t,x,A,H,N,B) A*x + H*(kron(x,x)) + (N*x + B)*u(t);
Jac_xdot = @(t,x,A,H,N,B) A + 2*H*(kron(x,speye(size(A,1)))) + N*u(t); 
optionsFOM = odeset('RelTol',1e-6,'AbsTol',1e-9,'Jac',Jac_xdot);
x0QB = 0.1*ones(nQB,1);
tic
[~, xQB] = ode15s(xdotQB, tSim, x0QB, optionsFOM, A, H, N, B);
simTimeFOMQB = toc;
yQB = C*xQB.';
subplot(2,3,i)
surf(x,tSim,xQB);
shading interp
xlabel('x (m)')
ylabel('time (s)')
zlabel('u(x,y)')
title (sprintf('Viscosity %g', nu))
i=i+1;
end

% Danish Rafiq
% This script solves the Quadratic-Bilinear version of 1D Burgers' equation
clear all; close all; clc;
%% simulation parameter
%Space
k = 500;             % No. of Spatial points
L = 5;               % Length of Spatial variable
x=linspace(0,L,k);   % Spatial vector 'x'
%Time
t0 = 0;
tEnd = 3;
dt = 0.01;
tSim = t0:dt:tEnd;

nu = 0.1750; % Viscosity

% excitation signals
u = @(t) 0.5*(cos(2*pi*t/10)+1);
% u = @(t) exp(-t);
% u = @(t) 1;
% u = @(t) 2*sin(2*pi*t/100);

%% Simulation of the QB-version of Burger model
[E,A,H,N,B,C,v0] = Burgers_Matrices_nonzero_ic(k,nu); % E=I
nQB = size(A,1);
xdotQB = @(t,x,A,H,N,B) A*x + H*(kron(x,x)) + (N*x + B)*u(t);
Jac_xdot = @(t,x,A,H,N,B) A + 2*H*(kron(x,speye(size(A,1)))) + N*u(t); 
optionsFOM = odeset('RelTol',1e-3,'AbsTol',1e-6,'Jac',Jac_xdot);
x0QB = zeros(length(A),1);
tic
[~, xQB] = ode15s(xdotQB, tSim, x0QB, optionsFOM, A, H, N, B);
simTimeFOMQB = toc
xtQB = xQB + ones(size(xQB,1),1)*v0.';
yQB = C*xtQB.';
figure;
surf(x,tSim,xtQB);
shading interp
set(gca,'FontSize',20,'TickLabelInterpreter','LaTeX')
xlabel('$x$ (m)', 'Interpreter', 'LaTeX')
ylabel('t (s)' ,'Interpreter', 'LaTeX')
zlabel('$u(x,t)$', 'Interpreter', 'LaTeX')
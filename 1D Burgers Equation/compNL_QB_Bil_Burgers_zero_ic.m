% Danish Rafiq
clear all; close all; clc;
%% simulation parameter
%Space
k = 1500;             % No. of Spatial points
L = 5;               % Length of Spatial variable
x=linspace(0,L,k);   % Spatial vector 'x'
nu=0.05;

% Time
t0 = 0;
tEnd = 8;
dt = 0.01;
tSim = t0:dt:tEnd;

% excitation signals
u = @(t) 0.5*(cos(2*pi*t/10)+1);
% u = @(t) exp(-t);
% u = @(t) 1;
% u = @(t) 2*sin(pi*t);
%% Simulation of the QB-version of Burger model
[E,A,H,N,B,C] = Burgers_Matrices_zero_ic(k,nu); % E=I

nQB = size(A,1);

% Function and Jacobian of rhs of SISO QB-FOM (E=I)
xdotQB = @(t,x,A,H,N,B) A*x + H*(kron(x,x)) + (N*x + B)*u(t);
Jac_xdot = @(t,x,A,H,N,B) A + 2*H*(kron(x,speye(size(A,1)))) + N*u(t); 
optionsFOM = odeset('RelTol',1e-6,'AbsTol',1e-9,'Jac',Jac_xdot);

x0QB = zeros(nQB,1);

tic
[~, xQB] = ode15s(xdotQB, tSim, x0QB, optionsFOM, A, H, N, B);
simTimeFOMQB = toc
yQB = C*xQB.';
figure(1);
surf(x,tSim,xQB);
shading interp
set(gca,'FontSize',20,'TickLabelInterpreter','LaTeX')
xlabel('$x$ (m)', 'Interpreter', 'LaTeX')
ylabel('$t$ (s)' ,'Interpreter', 'LaTeX')
zlabel('$u(x,t)$', 'Interpreter', 'LaTeX')

%Plot Output
figure(2);
plot(tSim, yQB);
set(gca,'FontSize',20,'TickLabelInterpreter','LaTeX')
title('Output','Interpreter','LaTeX')
xlabel('$x$ (m)', 'Interpreter', 'LaTeX')
ylabel('$y$', 'Interpreter', 'LaTeX')



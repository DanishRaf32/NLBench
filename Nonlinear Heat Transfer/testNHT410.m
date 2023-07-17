%% Benchmark collection of large-scale nonlinear state-space models for dimensionality reduction
% Authors: Danish Rafiq, M. A. Bazaz 
clear,clc,close all;

% System Matrices and Function Handle
% System: E x_dot + A x = B u + F f(x,u)
%         y = C x
[A, B, C, E, ~, f] = NHT410;
n = size(A,1);
%% Solving FOM
tic
dt = 1;
tEnd = 3000;
tSim = 0:dt:tEnd;

% excitation signal
u = @(t) [5e4; 0]; % u(t) = [J;Q] J: heat flux Q: heat source

x0 = zeros(n,1); % heat sink, T = 0 is room temperature

xdot = @(t,x) E\(B*u(t) + f(x) - A*x);

tic
[~,x] = ode23s(xdot,tSim,x0);
y = C*x.';
fulltime = toc;

figure;
plot(tSim, y);
xlabel('Time (s)');
ylabel('Outputs y');
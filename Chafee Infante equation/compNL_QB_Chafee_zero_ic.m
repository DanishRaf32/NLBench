clear all; close all; clc;

%% simulation parameter
k = 500;
L = 1;

t0 = 0;
tEnd = 5;
dt = 0.01;
tSim = t0:dt:tEnd;

% excitation signals
u = @(t) 0.5*(cos(pi*t)+1);
% u = @(t) 1;

%% Simulation of pure nonlinear Chafee-Infante model (combination: 2+3)
[E,A,B,C] = ChafeeInfante_POD_Matrices(k,L); %E=I

nNL = size(A,1)

f = fcubic3(k);

x0 = sparse(k,1);
xdot = @(t,x,A,B) A*x + B*u(t) - f(x);

Jac_xdot = @(t,x,A,B) A - ChafeeInfante_Jac2(x);
optionsFOM = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot);

tic
[~, xNL] = ode15s(xdot, tSim, x0, optionsFOM, A, B);
simTimeFOMNL = toc
yNL = C*xNL.';

figure;
plot(tSim,yNL);
xlabel('time ($s$)','Interpreter','LaTeX')
ylabel('Output $y(t)$','Interpreter','LaTeX')
grid on
set(gca,'FontSize',15,'TickLabelInterpreter','LaTeX')
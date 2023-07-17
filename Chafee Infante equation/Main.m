clear all
close all
clc 

%% Load benchmark
% Chafee Infante
% xdot = A*x + B*u(t) - fcubic(x)
k = 2000;
L = 1;
[E,A,B,C] = ChafeeInfante_POD_Matrices(k,L); %E=I

%% Simulations
% Definition of simulation time
t0 = 0;
tEnd = 2;
dt = 0.01;
tSim = t0:dt:tEnd;

% Definition of the simulation input
u = @(t) 0.5*(cos(2*pi*t)+1);
%u= @(t) 0.5*exp(-t)*(sin(pi*t)+1); 
%u= @(t) square(2*t);
%u = @(t) 1;

f = fcubic3(k);
fJac = ChafeeInfante_Jac3(k);
x0 = sparse(k,1);
xdot = @(t,x) A*x + B*u(t) - f(x);
Jac_xdot = @(t,x,A,B) A - fJac(x);
%Option for ode:
optionsFOM = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot);
%Option for Euler:
optionsNL.funJacobian = Jac_xdot;
optionsNL.Display = 'none';        % 'none', 'iter'
optionsNL.solver = 'fsolve';       % 'fsolve','NewtonRaphson'
optionsNL.RelTol = 1e-6;
optionsNL.AbsTol = 1e-9;
tic
[~, x] = implicitEuler(xdot, tSim, x0, optionsNL);
y = C*x.';
time.simulation.original = toc;

%% Output over space and time
discPoints = 0:L/k:(L-L/k);
figure(1)
surf(discPoints(1:5:end),tSim,x(:,1:5:end));
shading interp
grid on
xlabel('$\mathbf{x}$ (m)','Interpreter','LaTeX')
ylabel('$\mathbf{t} (s)$','Interpreter','LaTeX')
zlabel('$\vartheta(\mathbf{x},t)$','Interpreter','LaTeX')

%% Output over time
figure(2)
plot(tSim,y,'b');
grid on
legend('Full-order')
xlabel('Time(s)')
ylabel('v_1')

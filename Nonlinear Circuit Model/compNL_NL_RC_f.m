clear all; close all; clc;

%% simulation parameter
k = 500;
% k = 100000;
% k = 500000;

dt = 1e-2;
tEnd = 2;
tSim = 0:dt:tEnd;

% excitation signals
% u = @(t) 0*(t>=0);
% u = @(t) 0*(t<4) + 3*(t>=4);
% u = @(t) 0*(t<0.3) + 1*(t>=0.3);
% u = @(t) 2*exp(-t);
u = @(t) exp(-t);
% u = @(t) 30*exp(0.1*t);
% u = @(t) t.^2;
% u = @(t) 10*sin(3*exp(0.1*t)*t);
% u = @(t) 10*sin(t.^3);
% u = @(t) 10*3*exp(0.1*t);

% u = @(t) (cos(2*pi*t*10)-6)*2;
 
% f1 = 10; f2 = 1000;
% u = @(t) sin(2*pi*f1*t)+sin(2*pi*f2*t);
% u = @(t) 5*sin(2*pi*f1*t);
% u = @(t) sin(2*pi*f1*t).*sin(2*pi*f2*t);

%% Simulation of pure nonlinear RC-ladder model (with RC_ladder)
[f,B,C] = RC_ladder(k); % E=I
% xdot = RC_ladder(k); % E=I

B = zeros(k,1);
B(1) = 1;

C = zeros(1,k);
C(1) = 1;

nNL = k

xdotNL = @(t,x) f(x) + B*u(t);
% xdotNL = @(t,x) xdot(x,u(t));
optionsFOMNL = odeset('RelTol',1e-6,'AbsTol',1e-9,'Jacobian',@(t,x) RC_ladder_Jac(x,k));
% optionsFOMNL = odeset('RelTol',1e-6,'AbsTol',1e-9);

% x0NL = 0.1*ones(nNL,1);
x0NL = zeros(nNL,1);
tic
[~,xNL] = ode15s(xdotNL,tSim,x0NL,optionsFOMNL);
simTimeFOMNL = toc
yNL = C*xNL';

figure;
plot(tSim,yNL);
xlabel('Time (s)');
ylabel('Output y');

%% Simulation of pure nonlinear RC-ladder model (with RC_f)
xdotNL2 = @(t,x) RC_f(x) + B*u(t);

% optionsFOMNL = odeset('RelTol',1e-6,'AbsTol',1e-9);

% x0NL = 0.1*ones(nNL,1);
x0NL = zeros(nNL,1);
tic
[~,xNL2] = ode15s(xdotNL2,tSim,x0NL,optionsFOMNL);
simTimeFOMNL = toc
yNL2 = C*xNL2';

hold on;
plot(tSim,yNL2,'--');
legend('NL-FOM RC-ladder', 'NL-FOM RC-f');


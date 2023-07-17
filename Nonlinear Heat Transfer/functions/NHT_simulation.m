%% initialize and prepare
clear
clc
close all

%% system

%--------------------------------------------------------------------------
% parameters
%--------------------------------------------------------------------------
N=50;

a=[36 -0.1116 0.00017298 -1.78746e-7 1.3852815e-10];
L=0.1;
l=L/N;
A=1e-4;
rho=3970;
cp=766;

%--------------------------------------------------------------------------
% system
%--------------------------------------------------------------------------
% matrices
[E,K,B,C] = NHT_SystemMatrices(N,a,l,A,cp,rho);

% nonlinearity
f = @(x) NHT_nonlinearities_spd(x,a,N,A,l);

%% Simulation
%--------------------------------------------------------------------------
% simulation parameters
%--------------------------------------------------------------------------
dt = 1;
tEnd = 3000;
tSim = 0:dt:tEnd;

x0 = zeros(N,1);

% excitation signal
u = @(t) [5e4; 0]; % u(t) = [J;Q] J: heat flux; Q: heat source

%--------------------------------------------------------------------------
% definition of simulated system
%--------------------------------------------------------------------------
% [L_E,U_E] = lu(E);
% xdot = @(t,x) U_E\(L_E\(B*u(t) + f(x) - K*x));
% xdot = @(t,x) E\(B*u(t) + f(x) - K*x);
% [L,U,P,Q,R] = lu(E);            
% xdot = @(t,x) Q * (U \ (L \ (P * (R \ (B*u(t) + f(x) - K*x))))); 
xdot = @(t,x) (B*u(t) + f(x) - K*x);

%--------------------------------------------------------------------------
% simulation with or without analytically providing the Jacobian
%--------------------------------------------------------------------------
% Jac_xdot = @(t,x) U_E\(L_E\(NHT_Jac_spd(x,N,[],[],A,a,l) - K));
% Jac_xdot = @(t,x) E\(NHT_Jac_spd(x,N,[],[],A,a,l) - K);
% Jac_xdot = @(t,x) Q * (U \ (L \ (P * (R \ (NHT_Jac_spd(x,N,[],[],A,a,l) - K)))));

% Jac_xdot = @(t,x) NHT_Jac_spd(x,N,K,E,A,a,l);
Jac_xdot = @(t,x) (NHT_Jac_spd(x,N,[],[],A,a,l) - K);

optionsFOM = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot,'Mass',E);
% optionsFOM = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot);
% optionsFOM = odeset('Jacobian', Jac_xdot);
% optionsFOM.RelTol = 1e-6;

%--------------------------------------------------------------------------
% simulation
%--------------------------------------------------------------------------
tic
[~,x] = ode15s(xdot,tSim,x0,optionsFOM);
y = C*x';
simTime_yJac = toc

%--------------------------------------------------------------------------
% plotting the results
%--------------------------------------------------------------------------
% Output y
figure;
plot(tSim, y);

xlabel('Time [s]');
ylabel('Output y [K]');

% T(x,t=100)
X=0:1/N:1;
figure;
plot(X, [x(101,:) 0]);

xlabel('x/L [-]');
ylabel('T(x) [K]');
%% initialize and prepare
clear
clc
close all

%% system

%--------------------------------------------------------------------------
% parameters
%--------------------------------------------------------------------------
N=1000;

a1=[36 -0.1116 0.00017298 -1.78746e-7 1.3852815e-10];
a2=[36 -1.5 0.00017298 -1.78746e-7 1.3852815e-10];
a3=[36 -5 0.00017298 -1.78746e-7 1.3852815e-10];
L=0.1;
l=L/N;
A=1e-4;
rho=3970;
cp=766;

%% a1
[E,K,B,C] = NHT_SystemMatrices(N,a1,l,A,cp,rho);
f = @(x) NHT_nonlinearities_spd(x,a1,N,A,l);
%--------------------------------------------------------------------------
% Simulation
%--------------------------------------------------------------------------
% optionsJac = odeset('Jacobian',@(t,x) NHT_Jac_spd(x,N,K,E,A,a1,l),'Mass',E);
Jac_xdot = @(t,x) (NHT_Jac_spd(x,N,[],[],A,a1,l) - K);

optionsJac = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot,'Mass',E);
dt = 1;
tEnd = 3000;
tSim = 0:dt:tEnd;

% excitation signal
u = @(t) [5e4; 0]; % u(t) = [J;Q] J: heat flux Q: heat source

x0 = zeros(N,1);

xdot = @(t,x) (B*u(t) + f(x) - K*x);

tic
[~,x1] = ode15s(xdot,tSim,x0,optionsJac);
y1 = C*x1';
simTime_yJac = toc

%% a2
[E,K,B,C] = NHT_SystemMatrices(N,a2,l,A,cp,rho);
f = @(x) NHT_nonlinearities_spd(x,a2,N,A,l);
%--------------------------------------------------------------------------
% Simulation
%--------------------------------------------------------------------------
% optionsJac = odeset('Jacobian',@(t,x) NHT_Jac_spd(x,N,K,E,A,a2,l),'Mass',E);
Jac_xdot = @(t,x) (NHT_Jac_spd(x,N,[],[],A,a2,l) - K);
optionsJac = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot,'Mass',E);

xdot = @(t,x) (B*u(t) + f(x) - K*x);

tic
[~,x2] = ode15s(xdot,tSim,x0,optionsJac);
y2 = C*x2';
simTime_yJac = toc

%% a3
[E,K,B,C] = NHT_SystemMatrices(N,a3,l,A,cp,rho);
f = @(x) NHT_nonlinearities_spd(x,a3,N,A,l);
%--------------------------------------------------------------------------
% Simulation
%--------------------------------------------------------------------------
% optionsJac = odeset('Jacobian',@(t,x) NHT_Jac_spd(x,N,K,E,A,a3,l),'Mass',E);
Jac_xdot = @(t,x) (NHT_Jac_spd(x,N,[],[],A,a3,l) - K);
optionsJac = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot,'Mass',E);

xdot = @(t,x) (B*u(t) + f(x) - K*x);

tic
[~,x3] = ode15s(xdot,tSim,x0,optionsJac);
y3 = C*x3';
simTime_yJac = toc

%% Visualisation 

%-> T0 and T-middle
figure;
plot(tSim, y1(1,:), 'b');
hold on;
plot(tSim, y1(2,:), 'g');

plot(tSim, y2(1,:), 'r');
plot(tSim, y2(2,:), 'y');

plot(tSim, y3(1,:), 'm');
plot(tSim, y3(2,:), 'c');

xlabel('Time (s)');
ylabel('Output y [K]');

legend('a1, T(X=0)','a1, T(X=0.5L)','a2, T(X=0)','a2, T(X=0.5L)','a3, T(X=0)','a3, T(X=0.5L)')

%-> T(X) with X=normed location, t=const.
figure;
X=0:1/N:1;
plot(X, [x1(2,:) 0],'b -o');
hold on;
plot(X, [x1(50,:) 0],'b -+');
plot(X, [x1(300,:) 0],'b -*');
plot(X, [x1(1501,:) 0],'b -^');

plot(X, [x2(2,:) 0],'r -o');
plot(X, [x2(50,:) 0],'r -+');
plot(X, [x2(300,:) 0],'r -*');
plot(X, [x2(1501,:) 0],'r -^');

plot(X, [x3(2,:) 0],'g -o');
plot(X, [x3(50,:) 0],'g -+');
plot(X, [x3(300,:) 0],'g -*');
plot(X, [x3(1501,:) 0],'g -^');

xlabel('Location [-]');
ylabel('Temperature T(x) [K]');

legend('t=1','t=50','t=300','t=1500','t=1','t=50','t=300','t=1500','t=1','t=50','t=300','t=1500')


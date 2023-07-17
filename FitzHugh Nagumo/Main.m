clear all; close all; clc;

%% simulation parameter
k = 500;  % Order of system = 2*k
l=1;  % spatial length
t0 = 0;
tEnd = 15;
dt = 0.01;
tSim = t0:dt:tEnd;
% excitation signals
i0 = @(t) 5*10^4*t.^3.*exp(-15*t);
u2 = @(t) 1;
u = @(t) [i0(t); u2(t)];

%% Simulation of pure nonlinear FHN-model
[E,A,B,C] = Fitz_Matrices(k,l); % E~=I !!

A = E\A; %% Since E is diagonal, so these operations are very cheap and efficient. 
B = E\B; 
E = E\E; % E=I

Ed = E(1:2*k,1:2*k);
Ad = A(1:2*k,1:2*k);
Bd = B(1:2*k,:);
Cd = C(:,1:2*k);
nNL = size(A,1)
xdotNL = @(t,x,A,B) A*x + Fitz_NL(x) +  B*u(t);
options = odeset('RelTol',1e-10,'AbsTol',1e-12);
% options = odeset('RelTol',1e-10,'AbsTol',1e-12,'Mass',Ed);
x0NL = zeros(2*k,1);
tic
[~,xNL] = ode15s(xdotNL,tSim,x0NL,options,Ad,Bd);
simTimeFOMNL = toc
yNL = Cd*xNL';
%% Plotting
figure(1);
plot(tSim,yNL,'linewidth',2);
xlabel('time ($s$)','Interpreter','LaTex')
ylabel('Output','Interpreter','LaTex')
legend('y_{1}(t)','y_{2}(t)')
set(gca,'FontSize',15,'TickLabelInterpreter','LaTeX')
grid on
L=linspace(0,l,length(tSim));

figure(2)
plot3(L,yNL(1,:),yNL(2,:),'r','Linewidth',3)
grid on
title('Limit cycle behaviour','Interpreter','LaTeX')
set(gca,'FontSize',15,'TickLabelInterpreter','latex','XMinorGrid','on',...
    'YMinorGrid','on','ZMinorGrid','on');
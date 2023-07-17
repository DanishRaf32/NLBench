%% Benchmark collection of large-scale nonlinear state-space models for dimensionality reduction
% Authors: Danish Rafiq, M. A. Bazaz 

% Solves the Ring Grid Power System Model with N similar generators.
% Ref:

% [1] Structure preserving reduced order modeling technique for power
%     system (Rafiq D., and Bazaz. M. A., Seventh Indian Control Conference (ICC),2021)
% [2] Reduced order modelling for transient simulation of power systems
%     using TPWL(Haris Malik, Borzacchiello,Francisco and Diez)
clear; close all; clc
% System Parameters
N=1500;         % Full Order of the system
u_train=@(t) 0.52;
B=ones(2*N,1);
C=zeros(1,2*N);
C(1)=1;
IC=zeros(2*N,1); %initial conditions 
tsim=0:0.1:30;    % Simulation time

% Select Solver
solver_for_basis='ode15s';  %ode45s ode15s Euler
solver_for_ROM='ode15s';  %ode45s ode15s Euler

% Start Simulation of FOM with Selected Solver
f=@(x) Swing_equation(x);
xdotNL= @(t,x)f(x)+u_train(t);
Jack= @(t,x) SwingJack(x);
%Options for ode
optionsFOM = odeset('RelTol',1e-6,'AbsTol',1e-9,'Jac',Jack);
%Option for Euler
optionsNL.funJacobian = Jack;
optionsNL.Display = 'none';        % 'none', 'iter'
optionsNL.solver = 'fsolve';       % 'fsolve','NewtonRaphson'
optionsNL.RelTol = 1e-6;
optionsNL.AbsTol = 1e-9;
disp('Taking Snapshots of the system:')
switch(solver_for_basis)
    case 'ode45'
       tic 
       [tFOM,xFOM]=ode45(xdotNL,tsim,IC,optionsFOM);
       yFOM=C*xFOM';
       FOMSnaptime= toc
    case 'ode15s'
        tic
       [tFOM,xFOM]=ode15s(xdotNL,tsim,IC);
        yFOM=C*xFOM';
      FOMSnaptime=toc
    case 'Euler'
        tic
        [tFOM,xFOM,nNLSE] = implicitEuler(xdotNL,tsim,IC,optionsNL);
        yFOM=C*xFOM';
        FOMSnaptime=toc
    otherwise
        error('Please Select a solver')
end
plot(tsim,yFOM,'linewidth',3)
xlabel('t ($s$)','Interpreter','LaTeX')
ylabel('$\delta$','Interpreter','LaTeX')
grid on
hold on
set(gca,'FontSize',15,'TickLabelInterpreter','LaTeX')
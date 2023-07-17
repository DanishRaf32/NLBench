%% Benchmark collection of large-scale nonlinear state-space models for dimensionality reduction
% Authors: Danish Rafiq, M. A. Bazaz 
% File Crerated on 24.04.2021
%     Requirements:
      % 1. MATPOWER (https://matpower.org)
      % 2. pg_sync_models (https://sourceforge.net/projects/pg-sync-models/)

% Simulates IEEE 118 Power System Model via MATPOWER and pg_sync_models
%   [1] T. Nishikawa and A. E. Motter, Comparative analysis of existing
%       models for power-grid synchronization, New J. Phys. 17, 015012 (2015).
%
%   [2] A. E. Motter, S. A. Myers, M. Anghel, and T. Nishikawa, Spontaneous
%       synchrony in power-grid networks, Nat. Phys. 9, 191-197 (2013).
%
%% Initialize the IEEE Power System (Load data)
clear;clc;close all
mpc=case118;                  %118, 145, 300, 1888, 2736sp, test_system_10gen (see data file)
mpc.ref_freq=60;             % reference frequency 
[data,details]=EN_model(mpc); %EN, SM (H,D,A,K,gamma,omega_R)
n_oc=length(data.H);          %No of oscillators (FOM Size =2*n_oc)
tspan=0:0.01:8;

%% First Order System
f= @(x) power_func(x,data);
xdotNL= @(t,x) f(x);
IC=zeros(2*n_oc,1);
optionsNL.RelTol = 1e-6;
optionsNL.AbsTol = 1e-9;
%Option for Euler:
optionsNL.Display = 'none';        % 'none', 'iter'
optionsNL.solver = 'fsolve';       % 'fsolve','NewtonRaphson'
optionsNL.RelTol = 1e-6;
optionsNL.AbsTol = 1e-9;
disp('Simulating FOM...')
tic
[tFOM,xFOM]=implicitEuler(xdotNL,tspan,IC,optionsNL);
FOM_time=toc
C=zeros(1,2*n_oc);
C(n_oc+1:end)=1;
yFOM=C*xFOM'./n_oc;
delta=xFOM(:,1:n_oc);
FREQ = xFOM(:,n_oc+1:end);
%%
figure(1)
plot(tFOM,sum(delta,2)/n_oc,'linewidth',3);
grid on
ylabel('Avg. $ \delta$','Interpreter','LaTeX')
xlabel('$t (s)$','Interpreter','LaTeX')
set(gca,'FontSize',25,'TickLabelInterpreter','LaTeX')
set(gca,'FontSize',25,'TickLabelInterpreter','LaTeX')
set(gca,'FontSize',15,'TickLabelInterpreter','LaTeX')

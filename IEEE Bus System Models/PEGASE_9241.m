%% Benchmark collection of large-scale nonlinear state-space models for dimensionality reduction
% Authors: Danish Rafiq, M. A. Bazaz 
% File Crerated on 24.04.2021
%     Requirements:
      % 1. MATPOWER (https://matpower.org)
      % 2. pg_sync_models (https://sourceforge.net/projects/pg-sync-models/)

% Simulates Synthetic PEGASE 9241-bus power system model via MATPOWER and pg_sync_models
%     C. Josz, S. Fliscounakis, J. Maeght, and P. Panciatici, "AC Power Flow
%     Data in MATPOWER and QCQP Format: iTesla, RTE Snapshots, and PEGASE"
%     https://arxiv.org/abs/1603.01533
%
%     S. Fliscounakis, P. Panciatici, F. Capitanescu, and L. Wehenkel,
%     "Contingency ranking with respect to overloads in very large power
%     systems taking into account uncertainty, preventive and corrective
%     actions", Power Systems, IEEE Trans. on, (28)4:4909-4917, 2013.
%     https://doi.org/10.1109/TPWRS.2013.2251015
%
%    A. E. Motter, S. A. Myers, M. Anghel, and T. Nishikawa, Spontaneous
%       synchrony in power-grid networks, Nat. Phys. 9, 191-197 (2013).
%
%% Initialize the IEEE Power System (Load data)
clear;clc;close all
mpc=case9241pegase;             %Load data 
mpc.ref_freq=60;                %reference frequency
[data,details]=EN_model(mpc);   %EN, SM model 
n_oc=length(data.H);            %No of oscillators (FOM Size =2*n_oc)
tspan=0:0.01:10;

%% First Order System
f1= @(x) power_func(x,data);
xdotNL1= @(t,x) f1(x);
IC=zeros(2*n_oc,1);
disp('Simulating Power System Model...')
tic
[tFOM,xFOM]=ode15s(xdotNL1,tspan,IC);
FOM_time_EN=toc
C=zeros(1,2*n_oc);
C(n_oc+1:end)=1;
yFOM=C*xFOM'./n_oc;
delta=xFOM(:,1:n_oc);
FREQ = xFOM(:,n_oc+1:end);
%% Plotting
figure(1)
plot(tFOM,sum(delta,2)/n_oc,'linewidth',3);
grid on
ylabel('Avg. $\delta$','Interpreter','LaTeX')
xlabel('$t (s)$','Interpreter','LaTeX')
set(gca,'FontSize',20,'TickLabelInterpreter','LaTeX')
set(gca,'FontSize',20,'TickLabelInterpreter','LaTeX')
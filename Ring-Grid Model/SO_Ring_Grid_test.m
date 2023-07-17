%% Benchmark collection of large-scale nonlinear state-space models for dimensionality reduction
% Authors: Danish Rafiq, M. A. Bazaz  
% (created on 24.03.2021)
% Solves the Ring Grid Power System Model for 3 cases in second-order structure with n similar generators
% using implicit Euler method
% M*delta_ddot + D* delta_dot = pm-(b*sin(delta_i)+bint[sin(delta_i - delta_i+1)+sin(delta_i - delta_i-1)])
% Refs:
% [1] Structure preserving reduced order modeling technique for power
%     system (Rafiq D., and Bazaz. M. A., Seventh Indian Control Conference (ICC),2021)
% [2] Reduced order modelling for transient simulation of power systems
%     using TPWL(Haris Malik, Borzacchiello,Francisco and Diez)
clear; close all; clc
n=1500;         %Size of Full Order of the system
u=@(t)0.95;     %Pm=0.95
tf=35;          %Final time
del=1e-3;       %Step-size
tSim=0:del:tf-del;
k=round(tf/del);
%FOM 
M=eye(n);
D=0.25*eye(n);
B=ones(n,1);
C=ones(1,n);
f=@(x)-sin(x) - 100*(sin(x-circshift(x,1)) + sin(x-circshift(x,-1)));
Jf=@(x) -cos(x) - 100*(cos(x-circshift(x,1)) + cos(x-circshift(x,-1)));
%% Case-1: All nodes starting from under-perturbation
x0=0.8*ones(n,1); 
steps=k;
disp('Simulating Case-I: All nodes starting from under perturbation')
[xFOM1,FOM_time_case1]=implicitEulerSO(M,D,B,x0,f,del,u,steps);
yFOM1=C*xFOM1./n;  
FOM_time_case1
%% Case-II: All nodes starting from over-perturbation
x0=1.15*ones(n,1); 
steps=k;
disp('Simulating Case-II: All nodes starting from over-perturbation')
[xFOM2,FOM_time_case2]=implicitEulerSO(M,D,B,x0,f,del,u,steps);
yFOM2=C*xFOM2./n;  
FOM_time_case2
%% Case-III: All nodes starting from sychronously equillibrum condition
x0=ones(n,1); %(initial conditions) 0, 1, 1.15
steps=k;
disp('Simulating Case-III: All nodes starting from synchronously equillibrum condition')
[xFOM3,FOM_time_case3]=implicitEulerSO(M,D,B,x0,f,del,u,steps);
yFOM3=C*xFOM3./n; 
FOM_time_case3
%% Results 
disp('Case-III results:')
figure(1);
plot(tSim,yFOM1,'k',tSim,yFOM2,'r',tSim,yFOM3,'g','linewidth',3)
grid on
xlabel('t ($s$)','Interpreter','LaTeX')
ylabel('$\delta$','Interpreter','LaTeX')
legend('Case-I','Case-II','Case-III')
set(gca,'FontSize',15,'TickLabelInterpreter','LaTeX')


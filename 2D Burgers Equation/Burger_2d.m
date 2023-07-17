%% Benchmark collection of large-scale nonlinear state-space models for dimensionality reduction
% Author: Danish Rafiq, M. A. Bazaz
% Solves nonlinear viscous Burgers equation in 2D
% U_t + U*U_x + V*U_y = 1/Re (U_xx + U_yy)  [x,T] = [a,b]x[0,T]
% V_t + U*V_x + V*V_y = 1/Re (V_xx + V_yy)
% (x,y)=(a,b)x(c,d) , t=(0,T)
% subject to:
%Initial Cond:
% U(x,y,0) = Phi(x,y) and V(x,y,0) = phi(x,y) 
% and Boundary Conditions
% U(a,y,t)=g1; U(b,y,t)=g2; U(x,c,t)=g3; U(x,d,t)=g4;
% V(a,y,t)=g1; V(b,y,t)=g2; V(x,c,t)=g3; V(x,d,t)=g8;

clear all; clc; close all
a=0;b=1;c=0;d=1;
tend=1;
tsteps=250;
tspan=linspace(0,tend,tsteps);
Nx=30 ; Ny=30;             % Grid Size Total dimension: 2(Nx-2)(Ny-2)
Re=100;                    % Reynold' number
nxy=(Nx-2)*(Ny-2);
xx=linspace(a,b,Nx);
yy=linspace(c,d,Ny);
[Y,X]=meshgrid(xx,yy);
% Extract initial conditions from exact solution:
xx0=(3/4)-(1./(4*(1+exp((-4*X+4*Y)*Re/32))));
yy0=(3/4)+(1./(4*(1+exp((-4*X+4*Y)*Re/32))));
u0=xx0(2:Nx-1,2:Ny-1);
v0=yy0(2:Nx-1,2:Ny-1);
%Boundary condition for u and u form exact solution:
% For u
frame_u=zeros(Nx);
frame_u(:,1)=xx0(:,1);
frame_u(:,end)=xx0(:,end);
frame_u(1,2:end-1)=xx0(1,2:end-1);
frame_u(end,2:end-1)=xx0(end,2:end-1);
% For v
frame_v=zeros(Nx);
frame_v(:,1)=yy0(:,1);
frame_v(:,end)=yy0(:,end);
frame_v(1,2:end-1)=yy0(1,2:end-1);
frame_v(end,2:end-1)=yy0(end,2:end-1);
%%
x0=[reshape(u0',(Nx-2)*(Ny-2),1);reshape(v0',(Nx-2)*(Ny-2),1)];
fprintf('Simulating 2D Burger FD model...')
tic
[t,x]=ode15s(@(t,x)Burg_2d_FD(t,x,Nx,Ny,Re,a,b,c,d,xx0,yy0),tspan,x0); 
FOM_time=toc;
usol=x(:,1:nxy)';
vsol=x(:,nxy+1:end)';
% Plotting u at a given time step within tspan:
tt=10;    % t for plotting
%Save u
usol_tt=reshape(usol(:,tt),(Nx-2),(Nx-2))';
Usol_tt=frame_u;
Usol_tt(2:end-1,2:end-1)=usol_tt;
% Save v
vsol_tt=reshape(vsol(:,tt),(Nx-2),(Nx-2))';
Vsol_tt=frame_v;
Vsol_tt(2:end-1,2:end-1)=vsol_tt;
%Exact solution
tt1=tspan(tt);
Uexact_tt=(3/4)-(1./(4*(1+exp((-4*X+4*Y-tt1)*Re/32))));
Vexact_tt=(3/4)+(1./(4*(1+exp((-4*X+4*Y-tt1)*Re/32))));
%% Plot the numerical solution
figure(1)
set(gcf,'position',[20,20,1350,950])
subplot(2,2,1)
surf(xx,yy,Uexact_tt)
colormap('jet')
shading interp
colorbar
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex'); zlabel('$u(x,y)$','Interpreter','latex')
title(['Exact solution at time step $nt$=', num2str(tt)],'FontSize',25,'Interpreter','latex')
set(gca,'FontSize',25,'TickLabelInterpreter','latex','XMinorGrid','on',...
    'YMinorGrid','on');

subplot(2,2,2)
surf(Y,X,Usol_tt)
colormap('jet')
shading interp
colorbar
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex'); zlabel('$u(x,y)$','Interpreter','latex')
title(['Numerical solution at time step $nt$=', num2str(tt)],'FontSize',25,'Interpreter','latex')
set(gca,'FontSize',25,'TickLabelInterpreter','latex','XMinorGrid','on',...
    'YMinorGrid','on');

subplot(2,2,3)
surf(xx,yy,Vexact_tt)
colormap('jet')
shading interp
colorbar
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex'); zlabel('$v(x,y)$','Interpreter','latex')
title(['Exact solution at time step $nt$=', num2str(tt)],'FontSize',25,'Interpreter','latex')
set(gca,'FontSize',25,'TickLabelInterpreter','latex','XMinorGrid','on',...
    'YMinorGrid','on');

subplot(2,2,4)
surf(xx,yy,Vsol_tt)
colormap('jet')
shading interp
colorbar
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex'); zlabel('$v(x,y)$','Interpreter','latex')
title(['Numerical solution at time step $nt$=', num2str(tt)],'FontSize',25,'Interpreter','latex')
set(gca,'FontSize',25,'TickLabelInterpreter','latex','XMinorGrid','on',...
    'YMinorGrid','on');

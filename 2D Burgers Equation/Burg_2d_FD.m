function dxdt = Burg_2d_FD(t,x,Nx,Ny,Re,a,b,c,d,xx0,yy0)
deltax= (b-a)/(Nx-1); deltay= (d-c)/(Ny-1);
part=(Nx-2)*(Ny-2);
%part=numel(x)/2;
%xx=linspace(a,b,Nx);
%yy=linspace(c,d,Ny);
U=x(1:part);
V=x(part+1:end);
M1=sparse(zeros(Nx-2,Ny-2)+diag(ones((Nx-3),1),1)+diag(-1*ones((Nx-3),1),-1));
M=sparse(kron(eye(Nx-2),M1));
%N=(eye((Nx-2)*(Ny-2))-diag(ones((((Nx-2)*(Ny-2))-(Nx-2)),1),-(Nx-2)));
N=sparse(zeros((Nx-2)*(Ny-2))-diag(ones((((Nx-2)*(Ny-2))-(Nx-2)),1),-(Nx-2))+diag(ones((((Nx-2)*(Ny-2))-(Nx-2)),1),(Nx-2)));
f1=((1/(2*deltax))*M*U.*U) + ((1/(2*deltay))*N*U.*V);
f2=((1/(2*deltax))*M*V.*U) + ((1/(2*deltay))*N*V.*V);
D11=-2*eye(Nx-2)+diag(ones((Nx-3),1),1)+diag(ones((Nx-3),1),-1);
D1=sparse(kron(eye(Nx-2),D11));
D2=(-2*eye((Nx-2)*(Ny-2))+diag(ones((((Nx-2)*(Ny-2))-(Nx-2)),1),-(Nx-2))+diag(ones((((Nx-2)*(Ny-2))-(Nx-2)),1),(Nx-2)));
%Boundary Conditions
% For U: 
b1=zeros(Nx-2);
b1(:,1)=xx0(2:Nx-1,1);         % analytical
%b1(:,1)=ones(Nx-2,1)*0.5*exp(-t);  % external forcing
b1(:,end)=xx0(2:Nx-1,Ny);
   b1u=reshape(b1',(Nx-2)*(Ny-2),1);
b2=zeros(Nx-2);
b2(1,:)=xx0(1,2:Nx-1);
b2(end,:)=xx0(Ny,2:Nx-1);
   b2u=reshape(b2',(Nx-2)*(Ny-2),1);
B1=b1;
B1(:,end)=-B1(:,end);
B1=reshape(B1',(Nx-2)*(Ny-2),1);
   Bul=diag(B1);
B2=b2;
B2(end,:)=-B2(end,:);
B2=reshape(B2',(Nx-2)*(Ny-2),1);
   Bub=diag(B2);
% For V:
b11=zeros(Nx-2);
b11(:,1)=yy0(2:Nx-1,1);
%b11(:,1)=ones(Nx-2,1)*0.5*exp(-t);  % external forcing
b11(:,end)=yy0(2:Nx-1,Ny);
   b1v=reshape(b11',(Nx-2)*(Ny-2),1);
b22=zeros(Nx-2);
b22(1,:)=yy0(1,2:Nx-1);
b22(end,:)=yy0(Ny,2:Nx-1);
    b2v=reshape(b22',(Nx-2)*(Ny-2),1);
B11=b11;
B11(:,end)=-B11(:,end);
B11=reshape(B11',(Nx-2)*(Ny-2),1);
   Bvl=diag(B11);
B22=b22;
B22(end,:)=-B22(end,:);
B22=reshape(B22',(Nx-2)*(Ny-2),1);
   Bvb=diag(B22);
%% ode for U and V:
dxdt= [ ((1/(Re*deltax^2))*((D1*U)+b1u) + (1/(Re*deltay^2))*((D2*U)+b2u) + ((1/(2*deltax))*Bul*U) + ((1/(2*deltay))*Bub*V) -f1) ;...
    ((1/(Re*deltax^2))*((D1*V)+b1v) + (1/(Re*deltay^2))*((D2*V)+b2v) + ((1/(2*deltax))*Bvl*U) + ((1/(2*deltay))*Bvb*V) -f2) ];

end


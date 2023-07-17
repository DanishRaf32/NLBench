function [A,B,C,D,E,F,bb]=Burger_2d_FDM_matrices(xx0,yy0,Nx,Ny,a,b,c,d)
%nxy=(Nx-2)*(Ny-2);
deltax= (b-a)/(Nx-1); deltay= (d-c)/(Ny-1);
x=1/(2*deltax);
y=1/(2*deltay);
D11=-2*eye(Nx-2)+diag(ones((Nx-3),1),1)+diag(ones((Nx-3),1),-1);
D1=sparse(kron(eye(Nx-2),D11));
D2=(-2*eye((Nx-2)*(Ny-2))+diag(ones((((Nx-2)*(Ny-2))-(Nx-2)),1),-(Nx-2))+diag(ones((((Nx-2)*(Ny-2))-(Nx-2)),1),(Nx-2)));
M1=sparse(zeros(Nx-2,Ny-2)+diag(ones((Nx-3),1),1)+diag(-1*ones((Nx-3),1),-1));
M=sparse(kron(eye(Nx-2),M1));
N=sparse(zeros((Nx-2)*(Ny-2))-diag(ones((((Nx-2)*(Ny-2))-(Nx-2)),1),-(Nx-2))+diag(ones((((Nx-2)*(Ny-2))-(Nx-2)),1),(Nx-2)));
%Boundary Conditions
% For U: 
%b1u
b1=zeros(Nx-2);
%b1(:,1)=xx0(2:Nx-1,1);         % left boundary analytical
b1(:,1)=zeros(Nx-2,1);           %left boundary external forcing
b1(:,end)=xx0(2:Nx-1,Ny);       %right analytical
%b1(:,end)=zeros(Nx-2,1);        % right external
   b1u=reshape(b1',(Nx-2)*(Ny-2),1);
%bb1
bb1=zeros(Nx-2); % for no control
bb1(:,1)=ones(Nx-2,1); % for left boundary control
%bb1(:,end)=ones(Nx-2,1);   % for both sides
bb1=reshape(bb1',(Nx-2)*(Ny-2),1);

   
%b2u   
b2=zeros(Nx-2);
b2(1,:)=xx0(1,2:Nx-1);
b2(end,:)=xx0(Ny,2:Nx-1);
   b2u=reshape(b2',(Nx-2)*(Ny-2),1);
%Bul   
B1=b1;
B1(:,end)=-B1(:,end);
B1=reshape(B1',(Nx-2)*(Ny-2),1);
   Bul=diag(B1);
%Bub   
B2=b2;
B2(end,:)=-B2(end,:);
B2=reshape(B2',(Nx-2)*(Ny-2),1);
   Bub=diag(B2);
% For V:
%b1v
b11=zeros(Nx-2);
%b11(:,1)=yy0(2:Nx-1,1); % left analytical
b11(:,1)=zeros(Nx-2,1);  % left external forcing
b11(:,end)=yy0(2:Nx-1,Ny); %right analytical
%b11(:,end)=zeros(Nx-2,1); %right external
   b1v=reshape(b11',(Nx-2)*(Ny-2),1);
%for external forcing
bb2=zeros(Nx-2);          % for no boundary control
bb2(:,1)=ones(Nx-2,1);   % for left boundary control
%bb2(:,end)=ones(Nx-2,1); % both sides
bb2=reshape(bb2',(Nx-2)*(Ny-2),1);
%b2v
b22=zeros(Nx-2);
b22(1,:)=yy0(1,2:Nx-1);
b22(end,:)=yy0(Ny,2:Nx-1);
    b2v=reshape(b22',(Nx-2)*(Ny-2),1);
%Bvl    
B11=b11;
B11(:,end)=-B11(:,end);
B11=reshape(B11',(Nx-2)*(Ny-2),1);
   Bvl=diag(B11);
%Bvb
B22=b22;
B22(end,:)=-B22(end,:);
B22=reshape(B22',(Nx-2)*(Ny-2),1);
   Bvb=diag(B22);
 %A
 A=sparse(kron(eye(2),D1));
 %B
 B=sparse([b1u;b1v]);
 %C
 C=sparse(kron(eye(2),D2));
 %D
 D=sparse([b2u;b2v]);
 %E
 E=sparse([x*Bul y*Bub;x*Bvl y*Bvb]);
 %F
 F=sparse([x*M y*N; x*M y*N]);
 %bb
 bb=sparse([bb1;bb2]);
end
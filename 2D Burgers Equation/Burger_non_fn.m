function[f]= Burger_non_fn(x)
a=0;b=1;c=0;d=1; % Mesh domain
Nx=10 ; Ny=10;
deltax= (b-a)/(Nx-1); deltay= (d-c)/(Ny-1);
part=(Nx-2)*(Ny-2);
U=x(1:part);
V=x(part+1:end);
M1=sparse(zeros(Nx-2,Ny-2)+diag(ones((Nx-3),1),1)+diag(-1*ones((Nx-3),1),-1));
M=sparse(kron(eye(Nx-2),M1));
N=sparse(zeros((Nx-2)*(Ny-2))-diag(ones((((Nx-2)*(Ny-2))-(Nx-2)),1),-(Nx-2))+diag(ones((((Nx-2)*(Ny-2))-(Nx-2)),1),(Nx-2)));
f1=((1/(2*deltax))*M*U.*U) + ((1/(2*deltay))*N*U.*V);
f2=((1/(2*deltax))*M*V.*U) + ((1/(2*deltay))*N*V.*V);
f=[f1;f2];
end
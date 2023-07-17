function [ f,P ] = ESB_nonlinearities(n,s,Bndry,order,x,e0,l)
%calculation of the nonlinearities
% f= P*q, q= x(2n+1:4n,1);

P = sparse(4*n,4*n);
g = sparse(4*n,1);

k = 1/(4*pi*e0);

if Bndry==[1 n] 
    x = [0;x(1:2*n-2);0;x(2*n-1:4*n-2)]; %include zero nodes for Bndry=[1 n]
    Boundry=[1 2*n-1];
elseif Bndry==1
    x= [0;0;x(1:4*n-2)];
    Boundry=[1 2];
end

for i=1:1:2*n
    for j=1:1:2*n
        if i~=j
            % case (2)
            if(i<=n && j<=n) && not((i>=n+1 || sum(i==Bndry)>0)) && not((j>=n+1 || sum(j==Bndry)>0))
                P(2*n+i,2*n+j)=1/sqrt((i-j)^2*l^2+(x(2*i-1)-x(2*j-1))^2);
                 %Pdy (1)
                 g(2*i-1)=g(2*i-1)+(x(2*i-1)-x(2*j-1))/sqrt((i-j)^2*l^2+(x(2*i-1)-x(2*j-1))^2)^3*x(2*n+i)*x(2*n+j);
                 %Pdy (2)
                 g(2*j-1)=g(2*j-1)-(x(2*i-1)-x(2*j-1))/sqrt((i-j)^2*l^2+(x(2*i-1)-x(2*j-1))^2)^3*x(2*n+i)*x(2*n+j);
            end
            % case (3)
            if(i<=n && j<=n) && not(i>=n+1 || sum(i==Bndry)>0) && (j>=n+1 || sum(j==Bndry)>0)
                P(2*n+i,2*n+j)=1/sqrt((i-j)^2*l^2+x(2*i-1)^2);
                %Pdy (3)
                g(2*i-1)=g(2*i-1)+x(2*i-1)/sqrt((i-j)^2*l^2+x(2*i-1)^2)^3*x(2*n+i)*x(2*n+j);
            end
            % case (4)
            if(i<=n && j<=n) && (i>=n+1 || sum(i==Bndry)>0) && not(j>=n+1 || sum(j==Bndry)>0)
                P(2*n+i,2*n+j)=1/sqrt((i-j)^2*l^2+x(2*j-1)^2);
                %Pdy (4)
                g(2*j-1)=g(2*j-1)+x(2*j-1)/sqrt((i-j)^2*l^2+x(2*j-1)^2)^3*x(2*n+i)*x(2*n+j);
            end
            % case (8)
            if(i<=n && j>n) && not(i>=n+1 || sum(i==Bndry)>0) && i~=j-n
                P(2*n+i,2*n+j)=1/sqrt((j-n-i)^2*l^2+(x(2*i-1)+s)^2);
                %Pdy(5)
                g(2*i-1)=g(2*i-1)+(x(2*i-1)+s)/sqrt((j-n-i)^2*l^2+(x(2*i-1)+s)^2)^3*x(2*n+i)*x(2*n+j);
            end
            % case (9)
            if(j<=n && i>n) && not(j>=n+1 || sum(j==Bndry)>0) && i~=j+n
                P(2*n+i,2*n+j)=1/sqrt((i-n-j)^2*l^2+(x(2*j-1)+s)^2);
                %Pdy(6)
                g(2*j-1)=g(2*j-1)+(x(2*j-1)+s)/sqrt((i-n-j)^2*l^2+(x(2*j-1)+s)^2)^3*x(2*n+i)*x(2*n+j);
            end
            % case (10)
            if i<=n && not(i>=n+1 || sum(i==Bndry)>0) && j==i+n
                P(2*n+i,2*n+j)=1/abs(x(2*i-1)+s);
                %Pdy(7)
                g(2*i-1)=g(2*i-1)+sign(x(2*i-1)+s)/(x(2*i-1)+s)^2*x(2*n+i)*x(2*n+j);
            end
            % case (11)
            if j<=n && not(j>=n+1 || sum(j==Bndry)>0) && i==j+n
                P(2*n+i,2*n+j)=1/abs(x(2*j-1)+s);
                %Pdy(8)
               g(2*j-1)=g(2*j-1)+sign(x(2*j-1)+s)/(x(2*j-1)+s)^2*x(2*n+i)*x(2*n+j);
            end
        end
    end
end

P=k*P;
g=0.5*k*g;
f=-P*x+g;
%g=-Pdy*x; % -, f is moved from left side of = to the right side

% -------------------------------------------------------------------------
% Boundary Condition
% -------------------------------------------------------------------------

f(Boundry)=[];

%% first order system
if order == 1
    f = [sparse(length(f),1); f];
end
end
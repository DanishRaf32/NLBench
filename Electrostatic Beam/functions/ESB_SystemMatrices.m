function [varargout] = ESB_SystemMatrices(n,N,rho,E,I,l,h,w,A,s,e0,d1,d2,Bndry,order,Opts)
%calculates the system matrices: M ddz + D dz + K z = b*u + F*g
%n nodes are used per beam --> 2*N nodes
%N=n-1 elements per beam

% Syntax:
%       [M,D,K,b,c]  =  ESB_SystemMatrices(n,N,rho,E,I,l,h,w,A,s,e0,d1,d2,Bndry,2,'I')
%       [E,A,b,c]    =  ESB_SystemMatrices(n,N,rho,E,I,l,h,w,A,s,e0,d1,d2,Bndry,1,Opts)


Def.transf2nd = 'I'; %{ number, 'I', 'K', '-K' }

% create the options structure
if ~exist('Opts','var') || isempty(fieldnames(Opts))
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

%% calculation of the element mass- and stiffnes matrix
%--------------------------------------------------------------------------
% element stiffness matrix K
%--------------------------------------------------------------------------
Ke=sparse(4,4);
Ke=(2*E*I)/l^3 * [6 3*l -6 3*l; ...
                 3*l 2*l^2 -3*l l^2;...
                 -6 -3*l 6 -3*l;...
                 3*l l^2 -3*l 2*l^2];

%--------------------------------------------------------------------------
% element mass matrix
%--------------------------------------------------------------------------
Me_t=sparse(4,4);
Me_r=sparse(4,4);
Me_t=(rho*A*l)/420 * [156 22*l 54 -13*l;...
                    22*l 4*l^2 13*l -3*l^2;...
                    54 13*l 156 -22*l;...
                    -13*l -3*l^2 -22*l 4*l^2];
                
Me_r=(rho*I)/(30*l) * [36 3*l -36 3*l; ...
                     3*l 4*l^2 -3*l -l^2;...
                     -36 -3*l 36 -3*l; ...
                     3*l -l^2 -3*l 4*l^2];

Me=Me_t + Me_r;

%% assembly
K=sparse(4*n,4*n); %first node(1) at left pin and last node at right pin (N)
M=sparse(4*n,4*n);
for i = 1 : 1 : N
    a = 1 + 2*(i - 1) : 2*(i + 1);
    K(a,a) = K(a,a) + Ke;
    M(a,a) = M(a,a) + Me;
end

D = d1*M+d2*K;
%% constant entries of global capacity matrix (P) --> added to K
k=1/(4*pi*e0);
% -------------------------------------------------------------------------
% i==j --> diagonal entries
% -------------------------------------------------------------------------
K(2*n+1:4*n,2*n+1:4*n)=k*2/(w*l)*(l*log((w+sqrt(w^2+l^2))/l)+w*log((l+sqrt(w^2+l^2))/w))*speye(2*n,2*n);
% -------------------------------------------------------------------------
% all other entries need the value of i-j for calculation --> for ...
% -------------------------------------------------------------------------
for i=1:1:2*n
    for j=1:1:2*n
        if i~=j
            %5* (upper branch) || 5* (lower branch)
            if(i>n && j>n) || (i<=n && j<=n && (i>=n+1 || sum(i==Bndry)>0) && (j>=n+1||sum(j==Bndry)>0))     
                  K(2*n+i,2*n+j)=k*1/(abs(i-j)*l);
            end
            %6* (left branch)
            if (i<=n && j>n && (i>=n+1||sum(i==Bndry)>0))
                  K(2*n+i,2*n+j)=k*1/(sqrt((-i-n+j)^2*l^2+s^2));
            end
            %6* (right branch)
             if (i>n && j<=n && (j>=n+1||sum(j==Bndry)>0))
                  K(2*n+i,2*n+j)=k*1/(sqrt((i-n-j)^2*l^2+s^2));
             end
             %7* (left branch) || 7* (right branch)
            if (i<=n && (sum(i==Bndry)>0) && j==i+n) || (j<=n && (sum(j==Bndry)>0) && i==j+n)
                  K(2*n+i,2*n+j)=k*1/s;
            end
        end
    end    
end

% -------------------------------------------------------------------------
% b and c vector
% -------------------------------------------------------------------------
b = sparse(4*n,1);
b(2*n+1:3*n)=1;

c = sparse(1,4*n);

%% Boundary Conditions
if Bndry==[1 n]
    Boundry=[1 2*n-1];
    c(1,round(2*(n-1)/2)+1)=1;
elseif Bndry==1
    Boundry=[1 2];
    c(1,round(2*(n-1)/2)-1)=1;
end
K(:,Boundry) = [];
K(Boundry,:) = [];
D(:,Boundry) = [];
D(Boundry,:) = [];
M(:,Boundry) = [];
M(Boundry,:) = [];
b(Boundry) = [];
c(Boundry) = [];

%% first order system
% z = [x; xdot];
% E1st*zdot = A1st*z + b1st*u + f1st
if order == 1
    warning('The system is in 2nd order form and will be converted to 1st order.')
    
    % evaluate the specified option
    if isnumeric(Opts.transf2nd) 
        if isscalar(Opts.transf2nd)
            F = Opts.transf2nd * speye(size(K));
        else
            F = Opts.transf2nd;
        end 
    else 
        switch Opts.transf2nd 
            case 'I' 
                F = speye(size(K));
            case 'K' 
                F = K;
            case '-K' 
                F = -K;              
        end
    end
    
    E = [F sparse(size(K,1), size(K,1)); sparse(size(K,1), size(K,1)) M];
    A = [sparse(size(K,1), size(K,1)) F; -K -D];
    b = [sparse(length(b),1); b];
    c = [c sparse(1,length(c))];  
end

%% Output
if order == 1
    varargout{1}=E;
    varargout{2}=A;
    varargout{3}=b;
    varargout{4}=c;
else
    varargout{1}=M;
    varargout{2}=D;
    varargout{3}=K;
    varargout{4}=b;
    varargout{5}=c;
end
end



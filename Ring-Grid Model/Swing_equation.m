%Danish Rafiq
%Function for high-fidelity Swing Equation Model 
function [ f ] = Swing_equation(x)
bint=100;        % Susceptance between consecutive generators
d=0.45;          % Damping of the generators
%d=0.25;
b=1;             % Susceptance between generator and slack node
%pm=@(t) 0.95;    % Mechanical input
nnodes=numel(x)/2;      
f = zeros(size(x)); 
f(1:nnodes)=x(nnodes+1:end); %x1dot =x2;
delta = x(1:nnodes);   %x1
f(nnodes+1:end)=  - d*f(1:nnodes) - b*sin(delta) - ...
                bint*(sin(delta-circshift(delta,-1))+sin(delta-circshift(delta,1))); %x2dot           
end

%n=numel(x)/2;
%delta=x(1:n) %x1

%f=x(n+1:end) %x2=dotx1

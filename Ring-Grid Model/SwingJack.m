%Danish Rafiq
%Function for high-fidelity Swing Equation Model 
function [ Jf ] = SwingJack(x)
bint=100;        % Susceptance between consecutive generators
%d=0.25;          % Damping of the generators
b=1;             % Susceptance between generator and slack node
%pm=@(t) 0.95;    % Mechanical input
% nnodes=numel(x)/2;
% f = zeros(size(x));
% f(1:nnodes)=x(nnodes+1:end);
% delta = x(1:nnodes);  
% f(nnodes+1:end)=  -b*cos(delta)-bint*(cos(delta-circshift(delta,-1))+cos(delta-circshift(delta,1)));
 n=length(x);
 k=n/2;
 Jf = zeros(n);
% f(1:nnodes)=x(nnodes+1:end);
 delta = x(1:k);  
 Jf(1:k,1:k) =  diag(-b*cos(delta)-bint*(cos(delta-circshift(delta,-1))+cos(delta-circshift(delta,1))));
end
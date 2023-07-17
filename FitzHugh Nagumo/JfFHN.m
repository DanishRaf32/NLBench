function Jf = JacFitz_NL(x)
% Jacobian of nonlinear function Fitz_NL

n = length(x);
k = n/2; % n is always even
alpha = 1.1*(200/3);
beta = 200/3;

Jf = zeros(n);

for i = 1:k
    Jf(i,i) = 2*alpha*x(i) - 3*beta*x(i)^2;
end

end
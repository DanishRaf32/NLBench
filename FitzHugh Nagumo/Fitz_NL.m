function y = Fitz_NL(x)
% nonlinear function for FHN
n = length(x);
k = round(n/2);
alpha = 1.1*(200/3);
beta = 200/3;
y = alpha*[x(1:k).^2;sparse(k,1)] - beta*[x(1:k).^3;sparse(k,1)];
end
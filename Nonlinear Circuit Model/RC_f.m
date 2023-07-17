function y = RC_f(v)
% state-nonlinear function f for RC-Circuit
n = length(v);

y = zeros(n,1);

y(1) = - RC_g(v(1)) - RC_g(v(1)-v(2));

for i=2:n-1
   y(i) = RC_g(v(i-1)-v(i)) - RC_g(v(i)-v(i+1));    
end

y(n)= RC_g(v(n-1)-v(n));

end
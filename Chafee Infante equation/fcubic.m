function y = fcubic(v)

y=zeros(size(v));

for i = 1:length(v)
    y(i) = v(i)^3;
end

%%
% g = @(v) v.^3;
% k = length(v);
% A1 = spdiags(ones(k,k),0,k,k);
% y = g(A1*v);

end
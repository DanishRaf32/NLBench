function fJac = ChafeeInfante_Jac2(v)

gJac = @(v) 3*v.^2;
k = length(v);
A1 = speye(k);

fJac = transpose(A1)*spdiags(gJac(A1*v),0,k,k);

end
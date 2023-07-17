function fJac = ChafeeInfante_Jac3(k)

gJac = @(v) 3*v.^2;
A1 = speye(k);

fJac = @(v) transpose(A1)*spdiags(gJac(A1*v),0,k,k);

end
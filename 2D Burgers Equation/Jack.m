function fJac = Jack(nxy)

gJac = @(v) 2*v;
A1 = speye(2*nxy);

fJac = @(v) transpose(A1)*spdiags(gJac(A1*v),0,(2*nxy),(2*nxy));

end
% Unknowns elimination method (Master-slave)

function [jac, res] = ass_periodic_ms (strain_mac, u)

global dim
global X0Y0_nod
global X1Y0_nod
global X1Y1_nod
global X0Y1_nod

[jac, res] = ass_steff (u);

%corners
res([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) = 0.0; % x & y
res([X1Y0_nod*dim - 1, X1Y0_nod*dim + 0]) = 0.0; % x & y
res([X1Y1_nod*dim - 1, X1Y1_nod*dim + 0]) = 0.0; % x & y
res([X0Y1_nod*dim - 1, X0Y1_nod*dim + 0]) = 0.0; % x & y

jac([X0Y0_nod*dim - 1; X0Y0_nod*dim - 0], :) = 0.0;
jac([X1Y0_nod*dim - 1; X1Y0_nod*dim - 0], :) = 0.0;
jac([X1Y1_nod*dim - 1; X1Y1_nod*dim - 0], :) = 0.0;
jac([X0Y1_nod*dim - 1; X0Y1_nod*dim - 0], :) = 0.0;
jac(:,[X0Y0_nod*dim - 1; X0Y0_nod*dim - 0])  = 0.0;
jac(:,[X1Y0_nod*dim - 1; X1Y0_nod*dim - 0])  = 0.0;
jac(:,[X1Y1_nod*dim - 1; X1Y1_nod*dim - 0])  = 0.0;
jac(:,[X0Y1_nod*dim - 1; X0Y1_nod*dim - 0])  = 0.0;
for d = 0 : 1
 jac(X0Y0_nod*dim - d, X0Y0_nod*dim - d) = 1.0;
 jac(X1Y0_nod*dim - d, X1Y0_nod*dim - d) = 1.0;
 jac(X1Y1_nod*dim - d, X1Y1_nod*dim - d) = 1.0;
 jac(X0Y1_nod*dim - d, X0Y1_nod*dim - d) = 1.0;
end

endfunction

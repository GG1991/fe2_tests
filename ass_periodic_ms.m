
% Unknowns elimination method (Master-slave)

function [jac, res] = ass_periodic_ms (strain_mac, u_n)

global elements
global coordinates
global stress
global strain
global bc_nods
global bc_y0_per
global bc_y1_per
global bc_x0
global bc_x1
global lx
global ly
global dx
global dy
global nx
global ny
global nnods
global npe
global dim
global nelem
global size_tot

global X0Y0_nod
global X1Y0_nod
global X1Y1_nod
global X0Y1_nod

jac = sparse(size_tot, size_tot);
res = zeros(size_tot, 1);
u_e = zeros(npe*dim, 1);
ind = zeros(npe*dim, 1);

for e = 1 : nelem 
    u_e([1:2:npe*dim, 2:2:npe*dim]) = u_n([elements(e, :)*dim - 1, elements(e, :)*dim + 0]);

    [jac_e, res_e] = elemental (e, u_e);
    ind = [elements(e,:)*dim - 1; elements(e,:)*dim - 0](:);

    jac(ind, ind) += jac_e;
    res(ind) += res_e;
end

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

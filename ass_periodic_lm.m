% Lagrange multiplier method
function [jac, res] = ass_periodic_lm (strain_mac, u_n)

global elements
global coordinates
global stress
global strain
global bc_nods
global bc_y0
global bc_y1
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

jac = sparse(size_tot, size_tot);
res = zeros(size_tot, 1);
u_e = zeros(npe*dim, 1);
ind = zeros(npe*dim, 1);

[jac_1, res_1] = ass_steff (u_n);
jac([1:size(jac_1,1)],[1:size(jac_1,2)]) = jac_1;
res([1:size(res_1,1)]) = res_1;

X0Y0_nod = 1;
X1Y0_nod = nx;
X1Y1_nod = nx*ny;
X0Y1_nod = (ny-1)*nx + 1;

u_dif_y0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [0.0, ly]';
u_dif_x0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [lx , 0.0]';
u_X0Y0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [0.0, 0.0]';
u_X1Y0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [lx , 0.0]';
u_X1Y1 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [lx , ly ]';
u_X0Y1 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [0.0, ly ]';

res([1:2:max(size(bc_y0))*dim] + nx*ny*dim) = u_n([bc_y1*dim - 1]) - u_n([bc_y0*dim - 1]) - u_dif_y0(1); %x
res([2:2:max(size(bc_y0))*dim] + nx*ny*dim) = u_n([bc_y1*dim - 0]) - u_n([bc_y0*dim - 0]) - u_dif_y0(2); %y
res([1:2:max(size(bc_x0))*dim] + (nx*ny + max(size(bc_y0,2)) )*dim) = u_n([bc_x1*dim - 1]) - u_n([bc_x0*dim - 1]) - u_dif_x0(1); %x
res([2:2:max(size(bc_x0))*dim] + (nx*ny + max(size(bc_y0,2)) )*dim) = u_n([bc_x1*dim - 0]) - u_n([bc_x0*dim - 0]) - u_dif_x0(2); %y

% fb +- µ
res(bc_y0*dim - 1) -= u_n([1:2:max(size(bc_y0))*dim] + nx*ny*dim); %x
res(bc_y1*dim - 1) += u_n([1:2:max(size(bc_y0))*dim] + nx*ny*dim); %x
res(bc_y0*dim - 0) -= u_n([2:2:max(size(bc_y0))*dim] + nx*ny*dim); %y
res(bc_y1*dim - 0) += u_n([2:2:max(size(bc_y0))*dim] + nx*ny*dim); %y

res(bc_x0*dim - 1) -= u_n([1:2:max(size(bc_x0))*dim] + (nx*ny + max(size(bc_y0)))*dim); %x
res(bc_x1*dim - 1) += u_n([1:2:max(size(bc_x0))*dim] + (nx*ny + max(size(bc_y0)))*dim); %x
res(bc_x0*dim - 0) -= u_n([2:2:max(size(bc_x0))*dim] + (nx*ny + max(size(bc_y0)))*dim); %y
res(bc_x1*dim - 0) += u_n([2:2:max(size(bc_x0))*dim] + (nx*ny + max(size(bc_y0)))*dim); %y

res([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) = u_n([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) - u_X0Y0; % x & y
res([X1Y0_nod*dim - 1, X1Y0_nod*dim + 0]) = u_n([X1Y0_nod*dim - 1, X1Y0_nod*dim + 0]) - u_X1Y0; % x & y
res([X1Y1_nod*dim - 1, X1Y1_nod*dim + 0]) = u_n([X1Y1_nod*dim - 1, X1Y1_nod*dim + 0]) - u_X1Y1; % x & y
res([X0Y1_nod*dim - 1, X0Y1_nod*dim + 0]) = u_n([X0Y1_nod*dim - 1, X0Y1_nod*dim + 0]) - u_X0Y1; % x & y

for n = 1 : max(size(bc_y0))
  for d = 0 : 1
    jac(n*dim + nx*ny*dim - d, bc_y1(n)*dim - d) = +1.0;
    jac(n*dim + nx*ny*dim - d, bc_y0(n)*dim - d) = -1.0;
    jac(bc_y1(n)*dim - d, n*dim + nx*ny*dim - d) = +1.0;
    jac(bc_y0(n)*dim - d, n*dim + nx*ny*dim - d) = -1.0;
  end
end

for n = 1 : max(size(bc_x0))
  for d = 0 : 1
    jac(n*dim + ( nx*ny + max(size(bc_y0)) )*dim - d, bc_x1(n)*dim - d) = +1.0;
    jac(n*dim + ( nx*ny + max(size(bc_y0)) )*dim - d, bc_x0(n)*dim - d) = -1.0;
    jac(bc_x1(n)*dim - d, n*dim + (nx*ny + max(size(bc_y0)))*dim - d) = +1.0;
    jac(bc_x0(n)*dim - d, n*dim + (nx*ny + max(size(bc_y0)))*dim - d) = -1.0;
  end
end

jac([X0Y0_nod*dim - 1; X0Y0_nod*dim - 0], :) = 0.0;
jac([X1Y0_nod*dim - 1; X1Y0_nod*dim - 0], :) = 0.0;
jac([X1Y1_nod*dim - 1; X1Y1_nod*dim - 0], :) = 0.0;
jac([X0Y1_nod*dim - 1; X0Y1_nod*dim - 0], :) = 0.0;
for d = 0 : 1
 jac(X0Y0_nod*dim - d, X0Y0_nod*dim - d) = 1.0;
 jac(X1Y0_nod*dim - d, X1Y0_nod*dim - d) = 1.0;
 jac(X1Y1_nod*dim - d, X1Y1_nod*dim - d) = 1.0;
 jac(X0Y1_nod*dim - d, X0Y1_nod*dim - d) = 1.0;
end

endfunction

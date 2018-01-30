function [jac, res] = ass_periodic (strain_mac, u_n)

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


jac = sparse(size_tot, size_tot);
res = zeros(size_tot, 1);
u_e = zeros(npe*dim, 1);
ind = zeros(npe*dim, 1);

for e = 1 : nelem 

    u_e([1:2:npe*dim]) = u_n([elements(e, :)*dim - 1]); %set x vals
    u_e([2:2:npe*dim]) = u_n([elements(e, :)*dim + 0]); %set y vals

    [jac_e, res_e] = elemental (e, u_e);
    for n = 1 : npe 
      for d = 0 : 1
        ind(n*dim - d) = elements(e,n)*dim - d;
      end
    end

    jac(ind, ind) += jac_e;
    res(ind) += res_e;

end

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

res([1:2:max(size(bc_y0_per))*dim] + nx*ny*dim) = u_n([bc_y1_per*dim - 1]) - u_n([bc_y0_per*dim - 1]) - u_dif_y0(1); %x
res([2:2:max(size(bc_y0_per))*dim] + nx*ny*dim) = u_n([bc_y1_per*dim - 0]) - u_n([bc_y0_per*dim - 0]) - u_dif_y0(2); %y
res([1:2:max(size(bc_x0))*dim] + (nx*ny + max(size(bc_y0_per,2)) )*dim) = u_n([bc_x1*dim - 1]) - u_n([bc_x0*dim - 1]) - u_dif_x0(1); %x
res([2:2:max(size(bc_x0))*dim] + (nx*ny + max(size(bc_y0_per,2)) )*dim) = u_n([bc_x1*dim - 0]) - u_n([bc_x0*dim - 0]) - u_dif_x0(2); %y

% fb +- µ
res(bc_y0_per*dim - 1) -= u_n([1:2:max(size(bc_y0_per))*dim] + nx*ny*dim); %x
res(bc_y1_per*dim - 1) += u_n([1:2:max(size(bc_y0_per))*dim] + nx*ny*dim); %x
res(bc_y0_per*dim - 0) -= u_n([2:2:max(size(bc_y0_per))*dim] + nx*ny*dim); %y
res(bc_y1_per*dim - 0) += u_n([2:2:max(size(bc_y0_per))*dim] + nx*ny*dim); %y

res(bc_x0*dim - 1) -= u_n([1:2:max(size(bc_x0))*dim] + (nx*ny + max(size(bc_y0_per)))*dim); %x
res(bc_x1*dim - 1) += u_n([1:2:max(size(bc_x0))*dim] + (nx*ny + max(size(bc_y0_per)))*dim); %x
res(bc_x0*dim - 0) -= u_n([2:2:max(size(bc_x0))*dim] + (nx*ny + max(size(bc_y0_per)))*dim); %y
res(bc_x1*dim - 0) += u_n([2:2:max(size(bc_x0))*dim] + (nx*ny + max(size(bc_y0_per)))*dim); %y

res([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) = u_n([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) - u_X0Y0; % x & y
res([X1Y0_nod*dim - 1, X1Y0_nod*dim + 0]) = u_n([X1Y0_nod*dim - 1, X1Y0_nod*dim + 0]) - u_X1Y0; % x & y
res([X1Y1_nod*dim - 1, X1Y1_nod*dim + 0]) = u_n([X1Y1_nod*dim - 1, X1Y1_nod*dim + 0]) - u_X1Y1; % x & y
res([X0Y1_nod*dim - 1, X0Y1_nod*dim + 0]) = u_n([X0Y1_nod*dim - 1, X0Y1_nod*dim + 0]) - u_X0Y1; % x & y

for n = 1 : max(size(bc_y0_per))
  for d = 0 : 1
    jac(n*dim + nx*ny*dim - d, bc_y1_per(n)*dim - d) = +1.0;
    jac(n*dim + nx*ny*dim - d, bc_y0_per(n)*dim - d) = -1.0;
    jac(bc_y1_per(n)*dim - d, n*dim + nx*ny*dim - d) = +1.0;
    jac(bc_y0_per(n)*dim - d, n*dim + nx*ny*dim - d) = -1.0;
  end
end

for n = 1 : max(size(bc_x0))
  for d = 0 : 1
    jac(n*dim + ( nx*ny + max(size(bc_y0_per)) )*dim - d, bc_x1(n)*dim - d) = +1.0;
    jac(n*dim + ( nx*ny + max(size(bc_y0_per)) )*dim - d, bc_x0(n)*dim - d) = -1.0;
    jac(bc_x1(n)*dim - d, n*dim + (nx*ny + max(size(bc_y0_per)))*dim - d) = +1.0;
    jac(bc_x0(n)*dim - d, n*dim + (nx*ny + max(size(bc_y0_per)))*dim - d) = -1.0;
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

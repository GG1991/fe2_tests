function [jac, res] = ass_periodic (strain_mac, u_n)

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

u_dif_y0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [0.0, ly]';
u_dif_x0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [lx , 0.0]';

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

res([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) = u_n([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) - [0.0, 0.0]';

for n = 1 : size(bc_y0, 2)
  for d = 0 : 1
    jac(n*dim + nx*ny*dim - d, bc_y1(n)*dim - d) = +1.0;
    jac(n*dim + nx*ny*dim - d, bc_y0(n)*dim - d) = -1.0;
    jac(bc_y1(n)*dim - d, n*dim + nx*ny*dim - d) = +1.0;
    jac(bc_y0(n)*dim - d, n*dim + nx*ny*dim - d) = -1.0;
  end
end

for n = 1 : size(bc_x0, 2)
  for d = 0 : 1
    jac(n*dim + (nx*ny + size(bc_y0,2))*dim - d, bc_x1(n)*dim - d) = +1.0;
    jac(n*dim + (nx*ny + size(bc_y0,2))*dim - d, bc_x0(n)*dim - d) = -1.0;
    jac(bc_x1(n)*dim - d, n*dim + (nx*ny + size(bc_y0,2))*dim - d) = +1.0;
    jac(bc_x0(n)*dim - d, n*dim + (nx*ny + size(bc_y0,2))*dim - d) = -1.0;
  end
end

jac([X0Y0_nod*dim - 1; X0Y0_nod*dim - 0], :) = 0.0;
for d = 0 : 1
 jac(X0Y0_nod*dim - d, X0Y0_nod*dim - d) = 1.0;
end

endfunction

%  ensambla el Jacobiano (J) y el residuo (r)

function [jac, res] = ass_unifstrains (bc_nods, strain_mac, u_n)

global elements
global coordinates
global xg
global wg
global b_mat
global nx
global ny
global nnods
global lx
global ly
global dx
global dy
global npe
global dim
global nvoi

jac = sparse(dim*nnods, dim*nnods);
res = zeros(dim*nnods, 1);

for e = 1 : size(elements, 1) 

    u_e = u_n([elements(e, :)*dim - 1, elements(e, :)*dim]);
    [jac_e, res_e] = elemental (e, u_e);

    jac([elements(e,:)*dim - 1, elements(e,:)*dim], [elements(e,:)*dim - 1, elements(e,:)*dim]) += jac_e;
    res([elements(e,:)*dim - 1, elements(e,:)*dim]) += res_e;

end

u_d = zeros(size(bc_nods, 1)*dim, 1);
for n = 1 : size(bc_nods, 1)
  u_d([n*dim - 1, n*dim]) = [strain_mac(1) strain_mac(3) ; strain_mac(3) strain_mac(2)] * coordinates(bc_nods(n), :)';
end

res([bc_nods*dim - 1]) = u_n([bc_nods*dim]) - u_d([1:2:size(u_d,1)]); %set x vals
res([bc_nods*dim]) = u_n([bc_nods*dim]) - u_d([2:2:size(u_d,1)]);     %set y vals

jac([bc_nods*dim - 1; bc_nods*dim], :) = 0.0;
for n = 1 : size(bc_nods, 1)
  for d = 0 : 1
    jac(bc_nods(n)*dim - d, bc_nods(n)*dim - d) = 1.0;
  end
end

endfunction

%  ensambla el Jacobiano (J) y el residuo (r)

function [jac, res] = ass_unifstrains (strain_mac, u_n)

global elements
global coordinates
global bc_nods
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
global nelem

u_e = zeros(npe*dim,1);
jac = sparse(dim*nnods, dim*nnods);
res = zeros(dim*nnods, 1);

for e = 1 : nelem 

    u_e([1:2:npe*dim]) = u_n([elements(e, :)*dim - 1]); %set x vals
    u_e([2:2:npe*dim]) = u_n([elements(e, :)*dim + 0]); %set y vals

    [jac_e, res_e] = elemental (e, u_e);

    jac([elements(e,:)*dim - 1], [elements(e,:)*dim - 1]) += jac_e([1:2:size(jac_e,1)],[1:2:size(jac_e,1)]); %set x vals
    jac([elements(e,:)*dim + 0], [elements(e,:)*dim + 0]) += jac_e([2:2:size(jac_e,1)],[2:2:size(jac_e,1)]); %set y vals
    res([elements(e,:)*dim - 1]) += res_e([1:2:size(res_e,1)]); %set x vals
    res([elements(e,:)*dim + 0]) += res_e([2:2:size(res_e,1)]); %set y vals

end

u_d = zeros(size(bc_nods, 1)*dim, 1);
for n = 1 : size(bc_nods, 1)
  u_d([n*dim - 1, n*dim]) = [strain_mac(1) strain_mac(3) ; strain_mac(3) strain_mac(2)] * coordinates(bc_nods(n), :)';
end

res([bc_nods*dim - 1]) = u_n([bc_nods*dim - 1]) - u_d([1:2:size(u_d,1)]); %set x vals
res([bc_nods*dim + 0]) = u_n([bc_nods*dim + 0]) - u_d([2:2:size(u_d,1)]); %set y vals

jac([bc_nods*dim - 1; bc_nods*dim], :) = 0.0;
for n = 1 : size(bc_nods, 1)
  for d = 0 : 1
    jac(bc_nods(n)*dim - d, bc_nods(n)*dim - d) = 1.0;
  end
end

endfunction

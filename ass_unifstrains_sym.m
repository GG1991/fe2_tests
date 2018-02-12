function [jac, res] = ass_unifstrains_sym (strain_mac, u_n)

global elements
global coordinates
global bc_nods
global b_mat
global nnods
global npe
global dim
global nelem

u_e = zeros(npe*dim,1);
jac = sparse(dim*nnods, dim*nnods);
res = zeros(dim*nnods, 1);
ind = zeros(npe*dim, 1);

for e = 1 : nelem 
    u_e([1:2:npe*dim, 2:2:npe*dim]) = u_n([elements(e, :)*dim - 1, elements(e, :)*dim + 0]);
    [jac_e, res_e] = elemental (e, u_e);

    ind = [elements(e,:)*dim - 1; elements(e,:)*dim - 0](:);
    jac(ind, ind) += jac_e;
    res(ind) += res_e;
end

res([bc_nods*dim - 1, bc_nods*dim + 0]) = 0.0;

jac([bc_nods*dim - 1; bc_nods*dim], :) = 0.0;
jac(:, [bc_nods*dim - 1; bc_nods*dim]) = 0.0;
for n = 1 : size(bc_nods, 1)
  for d = 0 : 1
    jac(bc_nods(n)*dim - d, bc_nods(n)*dim - d) = 1.0;
  end
end

endfunction

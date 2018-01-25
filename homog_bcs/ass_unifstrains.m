function [jac, res] = ass_unifstrains (strain_mac, u_n)

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

u_d = zeros(size(bc_nods, 1)*dim, 1);
for n = 1 : size(bc_nods, 1)
  u_d([n*dim - 1, n*dim]) = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * coordinates(bc_nods(n), :)';
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

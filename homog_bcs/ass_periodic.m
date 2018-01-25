function [jac, res] = ass_periodic (strain_mac, u_n)

global elements
global coordinates
global bc_nods
global bc_y0
global bc_y1
global bc_x0
global bc_x1
global lx
global ly
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
    res([elements(e,:)*dim - 1]) += res_e([1:2:size(res_e,1)]); %set x vals
    res([elements(e,:)*dim + 0]) += res_e([2:2:size(res_e,1)]); %set y vals

end

u_dif_y0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [0.0, ly]';
u_dif_x0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [lx , 0.0]';

res([bc_y1*dim - 1]) = res([bc_y1*dim - 1]) + res([bc_y0*dim - 1]); %set x vals
res([bc_y1*dim - 0]) = res([bc_y1*dim - 0]) + res([bc_y0*dim - 0]); %set y vals
res([bc_x1*dim - 1]) = res([bc_x1*dim - 1]) + res([bc_x0*dim - 1]); %set x vals
res([bc_x1*dim - 0]) = res([bc_x1*dim - 0]) + res([bc_x0*dim - 0]); %set y vals

res([bc_y0*dim - 1]) = u_n([bc_y1*dim - 1]) - u_n([bc_y0*dim - 1]) - u_dif_y0(1); %set x vals
res([bc_y0*dim - 0]) = u_n([bc_y1*dim - 0]) - u_n([bc_y0*dim - 0]) - u_dif_y0(2); %set y vals
res([bc_x0*dim - 1]) = u_n([bc_x1*dim - 1]) - u_n([bc_x0*dim - 1]) - u_dif_x0(1); %set x vals
res([bc_x0*dim - 0]) = u_n([bc_x1*dim - 0]) - u_n([bc_x0*dim - 0]) - u_dif_x0(2); %set y vals

res([bc_y0(1)*dim - 1, bc_y0(1)*dim + 0]) = u_n([bc_y0(1)*dim - 1, bc_y0(1)*dim + 0]) - [0.0, 0.0]'

jac([bc_y0*dim - 1; bc_y0*dim - 0], :) = 0.0;
jac([bc_x0*dim - 1; bc_x0*dim - 0], :) = 0.0;
auxs = [0, 0];
for n = 1 : size(bc_y0, 2)
  for d = 0 : 1
    auxs = [jac(bc_y1(n)*dim - d, bc_y0(n)*dim - d), jac(bc_y1(n)*dim - d, bc_y1(n)*dim - d)];
    jac(bc_y0(n)*dim - d, bc_y1(n)*dim - d) = +1.0;
    jac(bc_y0(n)*dim - d, bc_y0(n)*dim - d) = -1.0;

    jac(bc_y1(n)*dim - d, :) = 0.0;

    jac(bc_y1(n)*dim - d, bc_y0(n)*dim - d) = auxs(1);
    jac(bc_y1(n)*dim - d, bc_y1(n)*dim - d) = auxs(2);
  end
end

for n = 1 : size(bc_x0, 2)
  for d = 0 : 1
    auxs = [jac(bc_x1(n)*dim - d, bc_x0(n)*dim - d), jac(bc_x1(n)*dim - d, bc_x1(n)*dim - d)];
    jac(bc_x0(n)*dim - d, bc_x1(n)*dim - d) = +1.0;
    jac(bc_x0(n)*dim - d, bc_x0(n)*dim - d) = -1.0;

    jac(bc_x1(n)*dim - d, :) = 0.0;

    jac(bc_x1(n)*dim - d, bc_x0(n)*dim - d) = auxs(1);
    jac(bc_x1(n)*dim - d, bc_x1(n)*dim - d) = auxs(2);
  end
end

jac([bc_y0(1)*dim - 1; bc_y0(1)*dim - 0], :) = 0.0;
for d = 0 : 1
 jac(bc_y0(1)*dim - d, bc_y0(1)*dim - d) = 1.0;
end

endfunction

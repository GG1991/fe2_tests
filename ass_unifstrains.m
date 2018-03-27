function [jac, res] = ass_unifstrains(strain_mac, u_n)

global bc_y0; global bc_y1; global bc_x0; global bc_x1; global dim
global X0Y0_nod; global X1Y0_nod; global X1Y1_nod; global X0Y1_nod

[jac, res] = ass_steff (u_n);

if(size(bc_x0,2) > 0)
 res([bc_x0*dim - 1; bc_x0*dim]) = 0.0;
 res([bc_x1*dim - 1; bc_x1*dim]) = 0.0;
endif

if(size(bc_y0,2) > 0)
 res([bc_y0*dim - 1; bc_y0*dim]) = 0.0;
 res([bc_y1*dim - 1; bc_y1*dim]) = 0.0;
endif

res([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) = 0.0; % x & y
res([X1Y0_nod*dim - 1, X1Y0_nod*dim + 0]) = 0.0; % x & y
res([X1Y1_nod*dim - 1, X1Y1_nod*dim + 0]) = 0.0; % x & y
res([X0Y1_nod*dim - 1, X0Y1_nod*dim + 0]) = 0.0; % x & y

if(size(bc_x0,2) > 0)
 jac([bc_x0*dim - 1; bc_x0*dim], :) = 0.0;
 jac([bc_x1*dim - 1; bc_x1*dim], :) = 0.0;
endif

if(size(bc_y0,2) > 0)
 jac([bc_y0*dim - 1; bc_y0*dim], :) = 0.0;
 jac([bc_y1*dim - 1; bc_y1*dim], :) = 0.0;
endif

if(size(bc_x0,2) > 0)
 jac(:, [bc_x0*dim - 1; bc_x0*dim]) = 0.0;
 jac(:, [bc_x1*dim - 1; bc_x1*dim]) = 0.0;
endif
if(size(bc_y0,2) > 0)
 jac(:, [bc_y0*dim - 1; bc_y0*dim]) = 0.0;
 jac(:, [bc_y1*dim - 1; bc_y1*dim]) = 0.0;
endif

jac([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0], :) = 0.0;
jac([X1Y0_nod*dim - 1, X1Y0_nod*dim + 0], :) = 0.0;
jac([X1Y1_nod*dim - 1, X1Y1_nod*dim + 0], :) = 0.0;
jac([X0Y1_nod*dim - 1, X0Y1_nod*dim + 0], :) = 0.0;
jac(:, [X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) = 0.0;
jac(:, [X1Y0_nod*dim - 1, X1Y0_nod*dim + 0]) = 0.0;
jac(:, [X1Y1_nod*dim - 1, X1Y1_nod*dim + 0]) = 0.0;
jac(:, [X0Y1_nod*dim - 1, X0Y1_nod*dim + 0]) = 0.0;

for d = 0 : 1

 if(size(bc_x0,2) > 0)
  for n = 1 : max(size(bc_x0))
    jac(bc_x0(n)*dim - d, bc_x0(n)*dim - d) = 1.0;
    jac(bc_x1(n)*dim - d, bc_x1(n)*dim - d) = 1.0;
  endfor
 endif

 if(size(bc_y0,2) > 0)
  for n = 1 : max(size(bc_y0))
    jac(bc_y0(n)*dim - d, bc_y0(n)*dim - d) = 1.0;
    jac(bc_y1(n)*dim - d, bc_y1(n)*dim - d) = 1.0;
  end
 endif

 jac(X0Y0_nod*dim - d, X0Y0_nod*dim - d) = 1.0;
 jac(X1Y0_nod*dim - d, X1Y0_nod*dim - d) = 1.0;
 jac(X1Y1_nod*dim - d, X1Y1_nod*dim - d) = 1.0;
 jac(X0Y1_nod*dim - d, X0Y1_nod*dim - d) = 1.0;
end

endfunction

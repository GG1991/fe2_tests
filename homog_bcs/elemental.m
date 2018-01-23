% fiber in the middle of a circle, structured mesh

function [jac_e, res_e] = elemental (e, elements, coordinates, lx, ly, xg, wg, dsh, b_mat, u_e)

dim = 2;
npe = 4;

jac_e = zeros(dim*npe, dim*npe);
res_e = zeros(dim*npe, 1);

nu = 0.3;
if ( distance(e, elements, coordinates, lx, ly) < 0.4 )
  E  = 1e7;
else
  E  = 1e6;
end
c_tan = [ 1-nu , nu   , 0           ;
          nu   , 1-nu , 0           ;
          0    , 0    , (1-2*nu)/2 ];

c_tan = c_tan * E / ((1 + nu)*(1 - 2*nu));

for gp = 1 : npe

  strain_gp = b_mat(:, :, gp) * u_e;
  stress_gp = c_tan * strain_gp;
  res_e += b_mat(:, :, gp)' * stress_gp;
  jac_e += b_mat(:, :, gp)' * c_tan * b_mat(:, :, gp);

end

endfunction

% fiber in the middle of a circle, structured mesh

function [jac_e, res_e] = elemental (e, u_e)

global elements
global coordinates
global elem_type
global lx
global ly
global xg
global wg
global b_mat
global dim
global npe
global strain
global stress
global mat_model

jac_e = zeros(dim*npe, dim*npe);
res_e = zeros(dim*npe, 1);


if (strcmp(mat_model,"plastic")) 

elseif (strcmp(mat_model,"linear")) 

  nu = 0.3;
  if ( distance(e) < 0.75 )
    elem_type(e) = 2;
    E  = 1e7;
  else
    elem_type(e) = 1;
    E  = 1e6;
  end
  c_tan = [ 1-nu , nu   , 0           ;
            nu   , 1-nu , 0           ;
            0    , 0    , (1-2*nu)/2 ];
  
  c_tan = c_tan * E / ((1 + nu)*(1 - 2*nu));
  
  strain(e,:) = [0.0 0.0 0.0];
  stress(e,:) = [0.0 0.0 0.0];
  
  for gp = 1 : npe
  
    strain_gp = b_mat(:, :, gp) * u_e;
    stress_gp = c_tan * strain_gp;
    strain(e,:) += strain_gp' * wg(gp);
    stress(e,:) += stress_gp' * wg(gp);
    res_e += b_mat(:, :, gp)' * stress_gp * wg(gp);
    jac_e += b_mat(:, :, gp)' * c_tan * b_mat(:, :, gp) * wg(gp);
  
  end

endif

strain(e,:) /= sum(wg);
stress(e,:) /= sum(wg);

endfunction

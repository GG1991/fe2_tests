% fiber in the middle of a circle, structured mesh

function [jac_e, res_e] = elemental (e, u_e)

global elements; global coordinates; global elem_type
global lx; global ly; global xg; global wg; global b_mat; global dim; global npe
global strain; global stress
global mat_model; global int_vars

jac_e = zeros(dim*npe, dim*npe);
res_e = zeros(dim*npe, 1);

nu = 0.3;
if ( distance(e) < 0.75 )
  elem_type(e) = 2;
  E  = 1e7;
else
  elem_type(e) = 1;
  E  = 1e6;
end

if (strcmp(mat_model,"plastic")) 
 
 E       = 1.0e9; 
 sig_y   = 2.0e11;
 d_eps   = 1.0e-6;

 for gp = 1 : npe

   eps_p_1 = int_vars((e-1)*4+gp, [1 2 3])';
   eps_e_1 = int_vars((e-1)*4+gp, [4 5 6])';
   eps_2 = b_mat(:, :, gp) * u_e;

   %calc c_tan by perturbations
   eps_2_1 = eps_2 + [d_eps; 0    ; 0    ];
   eps_2_2 = eps_2 + [0    ; d_eps; 0    ];
   eps_2_3 = eps_2 + [0    ; 0    ; d_eps];
   [sig_2_1, eps_e_dummy, eps_p_dummy] = model_plast(eps_2_1, eps_e_1, eps_p_1, E, nu, sig_y);
   [sig_2_2, eps_e_dummy, eps_p_dummy] = model_plast(eps_2_2, eps_e_1, eps_p_1, E, nu, sig_y);
   [sig_2_3, eps_e_dummy, eps_p_dummy] = model_plast(eps_2_3, eps_e_1, eps_p_1, E, nu, sig_y);
   c_tan(:,1) = sig_2_1/d_eps;
   c_tan(:,2) = sig_2_2/d_eps;
   c_tan(:,3) = sig_2_3/(2*d_eps);

   %calc sig_2
   [sig_2, eps_e_2, eps_p_2] = model_plast(eps_2, eps_e_1, eps_p_1, E, nu, sig_y);
   int_vars(e*4+(gp-1), [1 2 3]) = eps_p_2';
   int_vars(e*4+(gp-1), [4 5 6]) = eps_e_2';

   res_e += b_mat(:, :, gp)' * sig_2 * wg(gp);
   jac_e += b_mat(:, :, gp)' * c_tan * b_mat(:, :, gp) * wg(gp);

 endfor

elseif (strcmp(mat_model,"linear")) 

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
  endfor

endif

strain(e,:) /= sum(wg);
stress(e,:) /= sum(wg);

endfunction

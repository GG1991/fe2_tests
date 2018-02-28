% Unknowns elimination method (Master-slave)

function [jac_2, res_2] = ass_periodic_ms (strain_mac, u)

 global dim
 global X0Y0_nod
 global X1Y0_nod
 global X1Y1_nod
 global X0Y1_nod

 global ix_p
 global ix_m
 global ix_a
 
 [jac_1, res_1] = ass_steff (u);
 
 %corners
 res_1([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) = 0.0; % x & y
 res_1([X1Y0_nod*dim - 1, X1Y0_nod*dim + 0]) = 0.0; % x & y
 res_1([X1Y1_nod*dim - 1, X1Y1_nod*dim + 0]) = 0.0; % x & y
 res_1([X0Y1_nod*dim - 1, X0Y1_nod*dim + 0]) = 0.0; % x & y
 
 jac_1([X0Y0_nod*dim - 1; X0Y0_nod*dim - 0], :) = 0.0;
 jac_1([X1Y0_nod*dim - 1; X1Y0_nod*dim - 0], :) = 0.0;
 jac_1([X1Y1_nod*dim - 1; X1Y1_nod*dim - 0], :) = 0.0;
 jac_1([X0Y1_nod*dim - 1; X0Y1_nod*dim - 0], :) = 0.0;
 
 jac_1(:,[X0Y0_nod*dim - 1; X0Y0_nod*dim - 0])  = 0.0;
 jac_1(:,[X1Y0_nod*dim - 1; X1Y0_nod*dim - 0])  = 0.0;
 jac_1(:,[X1Y1_nod*dim - 1; X1Y1_nod*dim - 0])  = 0.0;
 jac_1(:,[X0Y1_nod*dim - 1; X0Y1_nod*dim - 0])  = 0.0;
 for d = 0 : 1
  jac_1(X0Y0_nod*dim - d, X0Y0_nod*dim - d) = 1.0;
  jac_1(X1Y0_nod*dim - d, X1Y0_nod*dim - d) = 1.0;
  jac_1(X1Y1_nod*dim - d, X1Y1_nod*dim - d) = 1.0;
  jac_1(X0Y1_nod*dim - d, X0Y1_nod*dim - d) = 1.0;
 end
  
 Kaa = jac_1(ix_a, ix_a); Kap = jac_1(ix_a, ix_p); Kam = jac_1(ix_a, ix_m); ra = res_1(ix_a);
 Kpa = jac_1(ix_p, ix_a); Kpp = jac_1(ix_p, ix_p); Kpm = jac_1(ix_p, ix_m); rp = res_1(ix_p); 
 Kma = jac_1(ix_m, ix_a); Kmp = jac_1(ix_m, ix_p); Kmm = jac_1(ix_m, ix_m); rm = res_1(ix_m); 
 jac_2 = [Kaa , (Kap+Kam); (Kma+Kpa), (Kpp+Kmp+Kpm+Kmm)];
 res_2 = [ra ; rm+rp];

endfunction

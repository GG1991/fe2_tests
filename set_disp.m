function [u] = set_disp(strain)

global dim; global bc_y0; global bc_y1; global bc_x0; global bc_x1;
global X0Y0_nod; global X1Y0_nod; global X1Y1_nod; global X0Y1_nod;
global bc_type; global size_tot; global coordinates; global lx; global ly;

u = zeros(size_tot, 1);

if (strcmp(bc_type,"ustrain"))

  if(size(bc_x0,2) > 0)
   for n = 1 : max(size(bc_x0))
     u([bc_x0(n)*dim - 1, bc_x0(n)*dim + 0]) = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * coordinates(bc_x0(n), :)';
     u([bc_x1(n)*dim - 1, bc_x1(n)*dim + 0]) = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * coordinates(bc_x1(n), :)';
   endfor
  endif
  if(size(bc_y0,2) > 0)
   for n = 1 : max(size(bc_y0))
     u([bc_y0(n)*dim - 1, bc_y0(n)*dim + 0]) = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * coordinates(bc_y0(n), :)';
     u([bc_y1(n)*dim - 1, bc_y1(n)*dim + 0]) = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * coordinates(bc_y1(n), :)';
   end
  endif
  u([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * [0.0, 0.0]'; % x & y
  u([X1Y0_nod*dim - 1, X1Y0_nod*dim + 0]) = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * [lx , 0.0]'; % x & y
  u([X1Y1_nod*dim - 1, X1Y1_nod*dim + 0]) = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * [lx , ly ]'; % x & y
  u([X0Y1_nod*dim - 1, X0Y1_nod*dim + 0]) = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * [0.0, ly ]'; % x & y

elseif (strcmp(bc_type,"per_ms"))

  u_X0Y0 = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * [0.0, 0.0]';
  u_X1Y0 = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * [lx , 0.0]';
  u_X1Y1 = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * [lx , ly ]';
  u_X0Y1 = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * [0.0, ly ]';
  
  % u+ = u- + c
  u_dif_y0 = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * [0.0, ly]';
  u_dif_x0 = [strain(1) strain(3)/2 ; strain(3)/2 strain(2)] * [lx , 0.0]';

  % u+ = u- + c
  u(bc_y1*dim - 1) = u(bc_y0*dim - 1) + u_dif_y0(1);
  u(bc_y1*dim - 0) = u(bc_y0*dim - 0) + u_dif_y0(2);
  u(bc_x1*dim - 1) = u(bc_x0*dim - 1) + u_dif_x0(1);
  u(bc_x1*dim - 0) = u(bc_x0*dim - 0) + u_dif_x0(2);
  
  % u_cor = eps * x
  u([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) = u_X0Y0; % x & y
  u([X1Y0_nod*dim - 1, X1Y0_nod*dim + 0]) = u_X1Y0; % x & y
  u([X1Y1_nod*dim - 1, X1Y1_nod*dim + 0]) = u_X1Y1; % x & y
  u([X0Y1_nod*dim - 1, X0Y1_nod*dim + 0]) = u_X0Y1; % x & y

endif

endfunction

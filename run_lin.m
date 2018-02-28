%
% solves micro problem for 1,2,3 experiments linear material
% can use bcs <ustrain|ustress|per_lm|per_ms>
% can use sol <lu|cg|cg_pd|cg_pgs>
%
global elements
global coordinates
global elem_type
global bc_nods
global bc_y0
global bc_y1
global bc_x0
global bc_x1
global xg
global wg
global b_mat
global stress
global strain
global res

global nx = 50;
global ny = 50;
global lx = 3;
global ly = 3;
global nn
global nelem
global dx
global dy

global npe = 4;
global dim = 2;
global nvoi = 3;

global X0Y0_nod
global X1Y0_nod
global X1Y1_nod
global X0Y1_nod

global ix_p
global ix_m
global ix_a

global size_tot;

init_vars();

bc_type   = "ustrain"; % <ustrain|ustress|per_lm|per_ms>
solver    = "lu";
for i = 1 : nargin
 if (strcmp(argv(){i}, "-cg"))
  solver = "cg";
 elseif (strcmp(argv(){i}, "-cg_pd"))
  solver = "cg_pd";
 elseif (strcmp(argv(){i}, "-cg_pgs"))
  solver = "cg_pgs";
 elseif (strcmp(argv(){i}, "-lu"))
  solver = "lu";
 elseif (strcmp(argv(){i}, "-ustrain"))
  bc_type = "ustrain";
  size_tot = nx*ny*dim;
 elseif (strcmp(argv(){i}, "-ustress"))
  bc_type = "ustress";
  size_tot = nx*ny*dim;
 elseif (strcmp(argv(){i}, "-per_lm"))
  bc_type = "per_lm";
  size_tot = (nx*ny + max(size(bc_y0)) + max(size(bc_x0))) * dim;
 elseif (strcmp(argv(){i}, "-per_ms"))
  bc_type = "per_ms";
  size_tot = nx*ny*dim;
 endif
end

printf("\033[33mBoundary Conditions = %s\n\033[0m",bc_type);
printf("\033[33mSolver = %s\n\033[0m",solver);

min_tol = 1.0e-7;
max_its = 3000;

c_ave = zeros(3,3);

strain_exp = [0.005 0 0; 0 0.005 0; 0 0 0.005]';

for i = 1 : 3

 printf ("\033[31mstrain = %f %f %f\n\033[0m", strain_exp(:,i)');
 
 u = zeros(size_tot, 1);

 if (bc_type == "per_ms")
   u_X0Y0 = [strain_exp(1,i) strain_exp(3,i)/2 ; strain_exp(3,i)/2 strain_exp(2,i)] * [0.0, 0.0]';
   u_X1Y0 = [strain_exp(1,i) strain_exp(3,i)/2 ; strain_exp(3,i)/2 strain_exp(2,i)] * [lx , 0.0]';
   u_X1Y1 = [strain_exp(1,i) strain_exp(3,i)/2 ; strain_exp(3,i)/2 strain_exp(2,i)] * [lx , ly ]';
   u_X0Y1 = [strain_exp(1,i) strain_exp(3,i)/2 ; strain_exp(3,i)/2 strain_exp(2,i)] * [0.0, ly ]';
   
   % u+ = u- + c
   u_dif_y0 = [strain_exp(1,i) strain_exp(3,i)/2 ; strain_exp(3,i)/2 strain_exp(2,i)] * [0.0, ly]';
   u_dif_x0 = [strain_exp(1,i) strain_exp(3,i)/2 ; strain_exp(3,i)/2 strain_exp(2,i)] * [lx , 0.0]';

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

 for nr = 1 : 3
 
   if (bc_type == "per_ms")
     [jac, res] = ass_periodic_ms (strain_exp(:,i), u);
   elseif (bc_type == "per_lm")
     [jac, res] = ass_periodic_lm (strain_exp(:,i), u);
   endif
 
   printf ("\033[32m|res| = %f\033[0m", norm(res));
   if (norm(res) < 1.0e-3)
     printf ("\n\033[0m");
     break 
   endif
 
   du = zeros(size(jac,2),1);

   tic()
    %[du, tol, its] = solver(jac, res, du, min_tol, max_its);
    if (strcmp(solver, "cg"))
     [du, tol, its] = cg(-jac, res, du, min_tol, max_its);
    elseif (strcmp(solver, "cg_pd"))
     [du, tol, its] = cg_pd(-jac, res, du, min_tol, max_its);
    elseif (strcmp(solver, "cg_pgs"))
     [du, tol, its] = cg_pgs(-jac, res, du, min_tol, max_its);
    elseif (strcmp(solver, "lu"))
     tol = 1; its = 1;
     du = -jac\res;
    endif
    
   if (bc_type == "per_ms")
    dua = du([1:size(ix_a,2)]);
    dum = du([size(ix_a,2) + 1 : size(du,1)]);
    dup = dum;
    u(ix_a) = u(ix_a) + dua;
    u(ix_m) = u(ix_m) + dum;
    u(ix_p) = u(ix_p) + dup;
   else
     u += du;
   endif

   time_sol = toc();
   printf ("\033[33m cg_tol = %f cg_its = %d cg_time = %f\n\033[0m", tol, its, time_sol);
 
 end

 [strain_ave, stress_ave] = average()
 c_ave(:,i) = stress_ave' / strain_ave(i);

end

printf ("\n");
c_ave

%figure();
%spy(jac); print -djpg spy.jpg 
%spy([Kaa , (Kap+Kam); (Kma+Kpa), (Kpp+Kmp+Kpm+Kmm)]); print -djpg spy.jpg 
%if ( issymmetric(full([Kaa , (Kap+Kam); (Kma+Kpa), (Kpp+Kmp+Kpm+Kmm)]),1.0e-8) )
%  printf ("\033[32mjac is symmetric\n");
%else
%  printf ("\033[31mjac is not symmetric\n");
%endif

[jac, res] = ass_periodic_ms (strain_exp(:,i), u);
write_vtk("sol.vtk", u)

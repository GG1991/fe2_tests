%
% solves micro problem for 1,2,3 experiments linear material
% 
% Usage:
% octave run_lin.m [solvers] [boundary conditions] [other running options]
% solvers               : [-lu|-cg|-cg_pd|-cg_pgs] 
% boundary conditions   : [-ustrain|-ustress|-per_ms|-per_lm]
% other running options : [-nexp 1]
%
global lx = 3; global ly = 3; global dx; global dy; global nn; global nelem; global nx; global ny;
global npe = 4; global dim = 2; global nvoi = 3; global size_tot;
global elements; global coordinates; global elem_type;global stress; global strain;
global bc_y0; global bc_y1; global bc_x0; global bc_x1;
global X0Y0_nod; global X1Y0_nod; global X1Y1_nod; global X0Y1_nod;
global ix_p; global ix_m; global ix_a; global xg; global wg; global b_mat;
global solver; global bc_type; global nexp;

% defaults
bc_type = "ustrain"; 
solver  = "lu";
nexp    = 3;
nx      = 5;
ny      = 5;

read_comm_line(nargin, argv);
init_vars();

printf("\033[33mBoundary Conditions = %s\n\033[0m",bc_type);
printf("\033[33mSolver = %s\n\033[0m",solver);
printf("\033[33mNumber of experiments = %d\n\033[0m",nexp);
printf("\033[33mUnknowns = %d\n\033[0m",size_tot);

min_tol = 1.0e-7;
max_its = 3000;

strain_exp = [0.005 0 0; 0 0.005 0; 0 0 0.005]';

for i = 1 : nexp

 printf ("\033[31mstrain = %f %f %f\n\033[0m", strain_exp(:,i)');
 
 u = set_disp(strain_exp(:,i));

 for nr = 1 : 3
 
   if (strcmp(bc_type,"ustrain"))
     [jac, res] = ass_unifstrains(strain_exp(:,i), u);
   elseif (strcmp(bc_type,"ustress"))
     [jac, res] = ass_unifstress_lm(strain_exp(:,i), u);
   elseif (strcmp(bc_type,"per_ms"))
     [jac, res] = ass_periodic_ms(strain_exp(:,i), u);
   elseif (strcmp(bc_type,"per_lm"))
     [jac, res] = ass_periodic_lm(strain_exp(:,i), u);
   endif
 
   printf ("\033[32m|res| = %f\033[0m", norm(res));
   if (norm(res) < 1.0e-3)
     printf ("\n\033[0m");
     break 
   endif
 
   du = zeros(size(jac,2),1);

   tic()
   if (strcmp(solver, "cg"))
    [du, tol, its] = cg(-jac, res, du, min_tol, max_its);
   elseif (strcmp(solver, "cg_pd"))
    [du, tol, its] = cg_pd(-jac, res, du, min_tol, max_its);
   elseif (strcmp(solver, "cg_pgs"))
    [du, tol, its] = cg_pgs(-jac, res, du, min_tol, max_its);
   elseif (strcmp(solver, "lu"))
    tol = 1; its = 1;
    du = -(jac\res);
   endif
   
   if (strcmp(bc_type,"per_ms"))
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
 if (nexp == 3) c_ave(:,i) = stress_ave' / strain_ave(i); end

end

if (nexp == 3) c_ave end

%figure();
%spy(jac); print -djpg spy.jpg 
%if ( issymmetric(full([Kaa , (Kap+Kam); (Kma+Kpa), (Kpp+Kmp+Kpm+Kmm)]),1.0e-8) )
%  printf ("\033[32mjac is symmetric\n");
%else
%  printf ("\033[31mjac is not symmetric\n");
%endif

[jac, res] = ass_periodic_ms (strain_exp(:,i), u);
write_vtk("sol.vtk", u)

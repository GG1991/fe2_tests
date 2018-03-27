%
% solves micro problem for 1,2,3 experiments linear material
% 
% Usage:
% octave run_lin.m [solvers] [boundary conditions] [other running options]
% solvers               : [-lu|-my_lu|-cg|-cg_pd|-cg_pgs] 
% boundary conditions   : [-ustrain|-ustress|-per_ms|-per_lm]
% other running options : [-nexp <n>]  sets number of experiments to <n>
%                         [-nx <n>]  <n> nodes in x direction
%                         [-ny <n>]  <n> nodes in y direction
%
global lx = 3; global ly = 3; global dx; global dy; global nn; global nelem; global nx; global ny;
global npe = 4; global dim = 2; global nvoi = 3; global size_tot;
global elements; global coordinates; global elem_type;global stress; global strain;
global bc_y0; global bc_y1; global bc_x0; global bc_x1;
global X0Y0_nod; global X1Y0_nod; global X1Y1_nod; global X0Y1_nod;
global ix_p; global ix_m; global ix_a; global xg; global wg; global b_mat;
global solver; global bc_type; global nexp;
global mat_model = 'plastic'; global int_vars;

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

dt = 0.01;
time_final = 0.35;
time_steps = round(time_final / dt);

strain_exp_0 = [0.005 0 0]';

for i = 1 : time_steps

 strain_exp = strain_exp_0 * (i-1) * dt;
 printf ("\033[31mstrain = %f %f %f\n\033[0m", strain_exp);
 u = set_disp(strain_exp);

 for nr = 1 : 3

   if (strcmp(bc_type,"ustrain"))
     [jac, res] = ass_unifstrains(strain_exp, u);
   elseif (strcmp(bc_type,"ustress"))
     [jac, res] = ass_unifstress_lm(strain_exp, u);
   elseif (strcmp(bc_type,"per_ms"))
     [jac, res] = ass_periodic_ms(strain_exp, u);
   elseif (strcmp(bc_type,"per_lm"))
     [jac, res] = ass_periodic_lm(strain_exp, u);
   end
 
   printf ("\033[32m|res| = %f\033[0m", norm(res));
   if (norm(res) < 1.0e-3) printf ("\n\033[0m"); break endif
 
   du = zeros(size(jac,2),1);

   tic()
   if (strcmp(solver, "cg"))
    [du, tol, its] = cg(-jac, res, du, min_tol, max_its);
   elseif (strcmp(solver, "cg_pd"))
    [du, tol, its] = cg_pd(-jac, res, du, min_tol, max_its);
   elseif (strcmp(solver, "cg_pgs"))
    [du, tol, its] = cg_pgs(-jac, res, du, min_tol, max_its);
   elseif (strcmp(solver, "cg_uzawa"))
   
     if (strcmp(bc_type,"per_lm"))
       A = jac([1:nx*ny*dim],[1:nx*ny*dim]);
       B = jac([1:nx*ny*dim],[nx*ny*dim+1 : (nx*ny + max(size(bc_y0)) + max(size(bc_x0)))*dim]);
       b1 = res([1:nx*ny*dim]);
       b2 = res([nx*ny*dim+1 : (nx*ny + max(size(bc_y0)) + max(size(bc_x0)))*dim]);
       x1 = zeros(nx*ny*dim, 1);
       x2 = zeros((max(size(bc_y0)) + max(size(bc_x0)))*dim, 1);
     end
     [x1, x2, tol, its] = cg_uzawa(A, B, b1, b2, x1, x2, min_tol, max_its);
     du = [x1;x2];

   elseif (strcmp(solver, "my_lu"))
    tol = 1; its = 1;
    du = my_lu(-jac, res);
   elseif (strcmp(solver, "lu"))
    tol = 1; its = 1;
    du = -(jac\res);
   end
   
   if (strcmp(bc_type,"per_ms"))
    dua = du([1:size(ix_a,2)]);
    dum = du([size(ix_a,2) + 1 : size(du,1)]);
    dup = dum;
    u(ix_a) = u(ix_a) + dua;
    u(ix_m) = u(ix_m) + dum;
    u(ix_p) = u(ix_p) + dup;
   else
    u += du;
   end

   time_sol = toc();
   printf ("\033[33m cg_tol = %f cg_its = %d cg_time = %f\n\033[0m", tol, its, time_sol);

 end

 %[strain_ave, stress_ave] = average()
 %if (nexp == 3) c_ave(:,i) = stress_ave' / strain_ave(i); end

end

%if (nexp == 3) c_ave end

%figure();
%spy(jac); print -djpg spy.jpg 
%if ( issymmetric(full([Kaa , (Kap+Kam); (Kma+Kpa), (Kpp+Kmp+Kpm+Kmm)]),1.0e-8) )
%  printf ("\033[32mjac is symmetric\n");
%else
%  printf ("\033[31mjac is not symmetric\n");
%endif
%write_vtk("sol.vtk", u)

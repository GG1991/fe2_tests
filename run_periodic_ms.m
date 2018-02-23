% Unknowns elimination method

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

global size_tot = nx*ny*dim;

init_vars();


strain_exp = [0.005 0 0; 0 0.005 0; 0 0 0.005]';

ix_p = [bc_x1*dim-1, bc_x1*dim-0, bc_y1*dim-1, bc_y1*dim-0]; % + indeces
ix_m = [bc_x0*dim-1, bc_x0*dim-0, bc_y0*dim-1, bc_y0*dim-0]; % + indeces
ix_a = setdiff([1:nn*dim],[ix_p, ix_m]); % interior indeces

c_ave = zeros(3,3);

for i = 1 : 3

u_X0Y0 = [strain_exp(1,i) strain_exp(3,i)/2 ; strain_exp(3,i)/2 strain_exp(2,i)] * [0.0, 0.0]';
u_X1Y0 = [strain_exp(1,i) strain_exp(3,i)/2 ; strain_exp(3,i)/2 strain_exp(2,i)] * [lx , 0.0]';
u_X1Y1 = [strain_exp(1,i) strain_exp(3,i)/2 ; strain_exp(3,i)/2 strain_exp(2,i)] * [lx , ly ]';
u_X0Y1 = [strain_exp(1,i) strain_exp(3,i)/2 ; strain_exp(3,i)/2 strain_exp(2,i)] * [0.0, ly ]';

% u+ = u- + c
u_dif_y0 = [strain_exp(1,i) strain_exp(3,i)/2 ; strain_exp(3,i)/2 strain_exp(2,i)] * [0.0, ly]';
u_dif_x0 = [strain_exp(1,i) strain_exp(3,i)/2 ; strain_exp(3,i)/2 strain_exp(2,i)] * [lx , 0.0]';

u = zeros(size_tot, 1);
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

printf ("\033[31mstrain = %f %f %f\n\033[0m", strain_exp(:,i)');

for nr = 1 : 3

  [jac, res] = ass_periodic_ms (strain_exp(:,i), u);

  Kaa = jac(ix_a, ix_a); Kap = jac(ix_a, ix_p); Kam = jac(ix_a, ix_m); ra = res(ix_a);
  Kpa = jac(ix_p, ix_a); Kpp = jac(ix_p, ix_p); Kpm = jac(ix_p, ix_m); rp = res(ix_p); 
  Kma = jac(ix_m, ix_a); Kmp = jac(ix_m, ix_p); Kmm = jac(ix_m, ix_m); rm = res(ix_m); 

  printf ("\033[32m|res| = %f\033[0m", norm([ra ; rm+rp]));
  if (norm([ra ; rm+rp]) < 1.0e-3)
    printf ("\n\033[0m");
    break 
  endif

  %du  = - [Kaa , (Kap+Kam); (Kma+Kpa), (Kpp+Kmp+Kpm+Kmm)] \ [ra ; rm+rp];
  du = zeros(size(Kaa,2) + size(Kap,2) ,1);
  min_tol = 1.0e-7;
  max_its = 3000;
  tic()
  tol = 1;
  its = 1;
% du = (-[Kaa , (Kap+Kam); (Kma+Kpa), (Kpp+Kmp+Kpm+Kmm)] ) \ [ra ; rm+rp]; tol = 0.0 ; its = 0;
% [du, tol, its]  = cg(-[Kaa , (Kap+Kam); (Kma+Kpa), (Kpp+Kmp+Kpm+Kmm)] , [ra ; rm+rp], du, min_tol, max_its);
  [du, tol, its]  = cg_pd(-[Kaa , (Kap+Kam); (Kma+Kpa), (Kpp+Kmp+Kpm+Kmm)] , [ra ; rm+rp], du, min_tol, max_its);
  time_sol = toc();
  printf ("\033[33m cg_tol = %f cg_its = %d cg_time = %f\n\033[0m", tol, its, time_sol);
  
  dua = du([1:size(ix_a,2)]);
  dum = du([size(ix_a,2) + 1 : size(du,1)]);
  dup = dum;

  u(ix_a) = u(ix_a) + dua;
  u(ix_m) = u(ix_m) + dum;
  u(ix_p) = u(ix_p) + dup;

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

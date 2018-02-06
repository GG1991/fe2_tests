
% Unknowns elimination method

global elements
global coordinates
global elem_type
global bc_nods
global bc_y0_per
global bc_y1_per
global bc_x0
global bc_x1
global xg
global wg
global b_mat
global stress
global strain
global res

global nx = 3;
global ny = 3;
global nelem = (nx-1)*(ny-1)
global nnods = nx*ny;
global lx = 3;
global ly = 3;
global dx = lx / (nx - 1);
global dy = ly / (ny - 1);

global npe = 4;
global dim = 2;
global nvoi = 3;

global size_tot = nx*ny*dim;

init_vars();

elem_type = zeros(nelem, 1);

du = zeros(size_tot, 1);
strain = zeros((nx-1)*(ny-1), nvoi);
stress = zeros((nx-1)*(ny-1), nvoi);

strain_exp = [0.005 0 0; 0 0.005 0; 0 0 0.005]';

c_ave = zeros(3,3);

for i = 1 : 1

u = zeros(size_tot, 1);
printf ("\033[31mstrain = %f %f %f\n\033[0m", strain_exp(:,i)');

for nr = 1 : 2

  [jac, res] = ass_periodic_ms (strain_exp(:,i), u);

  ix_p = [bc_x1*dim-1, bc_x1*dim-0, bc_y1_per*dim-1, bc_y1_per*dim-0]; % - indeces
  ix_m = [bc_x0*dim-1, bc_x0*dim-0, bc_y0_per*dim-1, bc_y0_per*dim-0]; % + indeces
  ix_a = setdiff([1:nx*ny*dim],[ix_p, ix_m]);
  
  Kaa = jac(ix_a, ix_a); % Kaa
  Kap = jac(ix_a, ix_p); % Ka+
  Kam = jac(ix_a, ix_m); % Ka-
  
  Kpa = Kap';
  Kpp = jac(ix_p, ix_p); 
  Kpm = jac(ix_p, ix_m);
  
  Kma = Kam';
  Kmp = Kpm'; 
  Kmm = jac(ix_m, ix_m);

  ra = res(ix_a);
  rp = res(ix_p);
  rm = res(ix_m);

  printf ("\033[32m|res| = %f\n\033[0m", norm(res));
  if (norm(res) < 1.0e-3); break; end

  du = - [Kaa (Kap + Kam) ; Kpa (Kpp + Kpm)] \ [ra - Kap*rm; rp - Kpp*rm ];
  dua = du([1:size(ix_a,2)]);
  dum = du([size(ix_a,2)+1:size(du,1)]);
  dup = rm + dum;

  dum 
  dup

  u(ix_a) = u(ix_a) + dua;
  u(ix_m) = u(ix_m) + dum;
  u(ix_p) = u(ix_p) + dup;

end

%[strain_ave, stress_ave] = average()
%c_ave(:,i) = stress_ave' / strain_ave(i);

end

%printf ("\n");
%c_ave

%figure();
%spy(jac); print -djpg spy.jpg 

write_vtk("sol.vtk", u)

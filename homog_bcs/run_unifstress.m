
% Lagrange multiplier Method (Non Linear Problem)

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

global nx = 40;
global ny = 40;
global nelem = (nx-1)*(ny-1)
global nnods = nx*ny;
global size_tot
global lx = 3;
global ly = 3;
global dx = lx / (nx - 1);
global dy = ly / (ny - 1);

global npe = 4;
global dim = 2;
global nvoi = 3;

init_vars();

size_tot = (nx*ny + nvoi) * dim;
elem_type = zeros(nelem, 1);

du = zeros(size_tot, 1);
strain = zeros((nx-1)*(ny-1), nvoi);
stress = zeros((nx-1)*(ny-1), nvoi);

strain_exp = [0.005 0 0; 0 0.005 0; 0 0 0.005]';

c_ave = zeros(3,3);

for i = 1 : 3

u = zeros(size_tot, 1);
printf ("\033[31mstrain = %f %f %f\n\033[0m", strain_exp(:,i)');

for nr = 1 : 2

  [jac, res] = ass_unifstress_lm (strain_exp(:,i), u);
  printf ("\033[32m|res| = %f\n\033[0m", norm(res));

  if (norm(res) < 1.0e-3)
    break
  end
  
  du = -(jac\res);
  u = u + du;

end

[strain_ave, stress_ave] = average()
c_ave(:,i) = stress_ave' / strain_ave(i);

end

printf ("\n");
c_ave

%figure();
%spy(jac); print -djpg spy.jpg 

write_vtk("sol.vtk", u)

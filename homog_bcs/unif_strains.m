global elements
global coordinates
global elem_type
global bc_nods
global xg
global wg
global b_mat
global stress
global strain

global nx = 4;
global ny = 4;
global nelem = (nx-1)*(ny-1)
global nnods = nx*ny;
global lx = 3;
global ly = 3;
global dx = lx / (nx - 1);
global dy = ly / (ny - 1);

global npe = 4;
global dim = 2;
global nvoi = 3;

init_vars();

elem_type = zeros(nelem, 1);

#elements
#coordinates
#bc_nods

u = zeros(nx*ny*dim, 1);
du = zeros(nx*ny*dim, 1);
strain = zeros((nx-1)*(ny-1), nvoi);
stress = zeros((nx-1)*(ny-1), nvoi);

dir_n = zeros(nx*ny*dim, 1);
strain_exp = [0.005 0 0; 0 0.005 0; 0 0 0.005]';

c_ave = zeros(3,3);

for i = 1 : 3

[jac, res] = ass_unifstrains (strain_exp(:,i), u);
printf ("\033[32m|res| = %f\n\033[0m", norm(res));

du = -(jac\res);
u = u + du;

[jac, res] = ass_unifstrains (strain_exp(:,i), u);
printf ("\033[32m|res| = %f\n\033[0m", norm(res));

[strain_ave, stress_ave] = average()
c_ave(:,i) = stress_ave' / strain_ave(i);

end

c_ave

%res

figure();
spy(jac); print -djpg spy.jpg 

write_vtk("sol.vtk", u)

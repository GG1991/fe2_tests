global elements
global coordinates
global elem_type
global bc_nods
global xg
global wg
global b_mat
global stress
global strain
global res

global nx = 50;
global ny = 50;
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

du = zeros(nnods*dim, 1);
strain = zeros((nx-1)*(ny-1), nvoi);
stress = zeros((nx-1)*(ny-1), nvoi);

dir_n = zeros(nx*ny*dim, 1);
strain_exp = [0.005 0 0; 0 0.005 0; 0 0 0.005];

c_ave = zeros(3,3);

tol = 0;
its = 0;

for i = 1 : 3

u = zeros(nnods*dim, 1);
u_d = zeros(size(bc_nods, 1)*dim, 1);
for n = 1 : size(bc_nods, 1)
  u_d([n*dim - 1, n*dim]) = [strain_exp(i,1) strain_exp(i,3)/2 ; strain_exp(i,3)/2 strain_exp(i,2)] * coordinates(bc_nods(n), :)';
end

u([bc_nods*dim - 1, bc_nods*dim + 0]) = u_d([1:2:size(u_d,1), 2:2:size(u_d,1)]);
printf ("\033[31mstrain = %f %f %f\n\033[0m", strain_exp(:,i)');

for nr = 1 : 2
  [jac, res] = ass_unifstrains_sym (strain_exp(:,i), u);
  printf ("\033[32m|res| = %f\n\033[0m", norm(res));
  if (norm(res) < 1.0e-3); break; end

  tic();
  %du = -(jac\res);
  [du, tol, its] = cg(jac, -res, du);
  time = toc()*1000;
  u = u + du;
end
printf("tol  :%e\n",tol);
printf("its  :%d\n",its);
printf("time :%f ms\n",time);

[strain_ave, stress_ave] = average()
c_ave(:,i) = stress_ave' / strain_ave(i);

end

printf ("\n");
c_ave

%figure();
%spy(jac); print -djpg spy.jpg 

write_vtk("sol.vtk", u)

global elements
global coordinates
global xg
global wg
global b_mat

global nx = 3;
global ny = 3;
global nnods = nx*ny;
global lx = 3;
global ly = 3;
global dx = lx / (nx - 1);
global dy = ly / (ny - 1);

global npe = 4;
global dim = 2;
global nvoi = 3;

init_vars();

#elements
#coordinates

u_n = zeros(nx*ny*dim, 1);
u = zeros(nx*ny*dim, 1);
du = zeros(nx*ny*dim, 1);

bc_y0 = [1 : 1 : nx];
bc_y1 = [(ny-1)*nx + 1 : 1 : nx*ny];
bc_x0 = [nx + 1 : nx : (ny-2)*nx + 1];
bc_x1 = [2*nx : nx : (ny-1)*nx];
bc_nods = [bc_y0, bc_y1, bc_x0, bc_x1]';

dir_n = zeros(nx*ny*dim, 1);
strain_exp = [0.005 0 0; 0 0.005 0; 0 0 0.005]';

[jac, res] = ass_unifstrains (bc_nods, strain_exp(:,1), u);
printf ("\033[32m|res| = %f\n\033[0m", norm(res));
res

du = jac\(-res);
u = u_n + du;

[jac, res] = ass_unifstrains (bc_nods, strain_exp(:,1), u);
printf ("\033[32m|res| = %f\n\033[0m", norm(res));
res

figure();
spy(jac); print -djpg spy.jpg 
quiver(coordinates(:,1), coordinates(:,2), u(1:2:nx*ny*2), u(2:2:nx*ny*2)); print -djpg sol.jpg


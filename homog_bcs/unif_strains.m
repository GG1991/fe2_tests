global elements
global coordinates
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

#elements
#coordinates
#bc_nods

u = zeros(nx*ny*dim, 1);
du = zeros(nx*ny*dim, 1);
strain = zeros((nx-1)*(ny-1), nvoi);
stress = zeros((nx-1)*(ny-1), nvoi);

dir_n = zeros(nx*ny*dim, 1);
strain_exp = [0.005 0 0; 0 0.005 0; 0 0 0.005]';

[jac, res] = ass_unifstrains (strain_exp(:,1), u);
printf ("\033[32m|res| = %f\n\033[0m", norm(res));

du = -(jac\res);
u = u + du;

[jac, res] = ass_unifstrains (strain_exp(:,1), u);
printf ("\033[32m|res| = %f\n\033[0m", norm(res));

du = -(jac\res);
u = u + du;

[jac, res] = ass_unifstrains (strain_exp(:,1), u);
printf ("\033[32m|res| = %f\n\033[0m", norm(res));

du = -(jac\res);
u = u + du;

[jac, res] = ass_unifstrains (strain_exp(:,1), u);
printf ("\033[32m|res| = %f\n\033[0m", norm(res));

du = -(jac\res);
u = u + du;

[jac, res] = ass_unifstrains (strain_exp(:,1), u);
printf ("\033[32m|res| = %f\n\033[0m", norm(res));

du = -(jac\res);
u = u + du;

[jac, res] = ass_unifstrains (strain_exp(:,1), u);
printf ("\033[32m|res| = %f\n\033[0m", norm(res));
du = -(jac\res);
u = u + du;

[jac, res] = ass_unifstrains (strain_exp(:,1), u);
printf ("\033[32m|res| = %f\n\033[0m", norm(res));

du = -(jac\res);
u = u + du;

[jac, res] = ass_unifstrains (strain_exp(:,1), u);
printf ("\033[32m|res| = %f\n\033[0m", norm(res));

[strain_ave, stress_ave] = average();

%res

%strain
%stress

%figure();
%spy(jac); print -djpg spy.jpg 
%quiver(coordinates(:,1), coordinates(:,2), u(1:2:nx*ny*2), u(2:2:nx*ny*2)); print -djpg sol.jpg

%x = coordinates(1:(nx-1)*(ny-1),1);
%y = coordinates(1:(nx-1)*(ny-1),2);
%z = stress(:,1)
%size(x)
%size(y)
%size(z)
%n = 9;
%[X, Y] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
%Z = griddata(x,y,z,X,Y);
%m = min(Z(Z~=0));
%M = max(Z(Z~=0));
%imshow((Z-m)/(M-m)); print -djpg map.jpg
%imshow(Z); print -djpg map.jpg

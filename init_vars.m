function init_vars()

global bc_nods; global bc_y0; global bc_y1; global bc_x0; global bc_x1
global xg; global wg; global b_mat; global dsh
global npe; global dim; global nvoi
global nx; global ny; global nn; global lx; global ly; global dx; global dy; global nelem
global X0Y0_nod; global X1Y0_nod; global X1Y1_nod; global X0Y1_nod
global ix_p; global ix_m; global ix_a; global solver; global bc_type; global size_tot
global elem_type; global strain; global stress; global elements; global int_vars;

nn = nx*ny;
nelem = (nx-1)*(ny-1);
dx = lx/(nx-1);
dy = ly/(ny-1);

elem_type = zeros(nelem, 1);
strain = zeros((nx-1)*(ny-1), nvoi);
stress = zeros((nx-1)*(ny-1), nvoi);
elements = zeros((nx-1)*(ny-1), npe);
int_vars = zeros(nelem*4, 7); % int_vars(e, [ eps_p_1(3) eps_e_1(3) alpha_1 ])

for i = 1 : (ny-1)
  for j = 1 : (nx-1)
    elements((i-1)*(nx-1) + j, 1) = j    + (i-1)*nx;
    elements((i-1)*(nx-1) + j, 2) = j+1  + (i-1)*nx;
    elements((i-1)*(nx-1) + j, 3) = j+1  + i*nx;
    elements((i-1)*(nx-1) + j, 4) = j    + i*nx;
  end
end

X0Y0_nod = 1;
X1Y0_nod = nx;
X1Y1_nod = nx*ny;
X0Y1_nod = (ny-1)*nx + 1;

bc_y0 = [2 : 1 : nx-1];
bc_y1 = [(ny-1)*nx + 2 : 1 : nx*ny-1];
bc_x0 = [nx + 1 : nx : (ny-2)*nx + 1];
bc_x1 = [2*nx : nx : (ny-1)*nx];
bc_nods = [bc_y0, bc_y1, bc_x0, bc_x1]';

if (strcmp(bc_type, "ustrain"))
 size_tot = nx*ny*dim;
elseif (strcmp(bc_type, "ustress"))
 size_tot = nx*ny*dim + nvoi;
elseif (strcmp(bc_type, "per_lm"))
 size_tot = (nx*ny + max(size(bc_y0)) + max(size(bc_x0))) * dim;
elseif (strcmp(bc_type, "per_ms"))
 size_tot = nx*ny*dim;
end

ix_p = [bc_x1*dim-1, bc_x1*dim-0, bc_y1*dim-1, bc_y1*dim-0]; % + indeces
ix_m = [bc_x0*dim-1, bc_x0*dim-0, bc_y0*dim-1, bc_y0*dim-0]; % + indeces
ix_a = setdiff([1:nn*dim],[ix_p, ix_m]); % interior indeces

global coordinates = zeros(nx*ny, dim);
for i = 1 : ny
  for j = 1 : nx
    coordinates((i - 1)*nx + j, 1) = (j - 1)*dx;
    coordinates((i - 1)*nx + j, 2) = (i - 1)*dy;
  end
end

wg = [0.25 , 0.25, 0.25 , 0.25]*(dx*dy);

xg = [ -0.577350269189626, -0.577350269189626;
       +0.577350269189626, -0.577350269189626;
       +0.577350269189626, +0.577350269189626;
       -0.577350269189626, +0.577350269189626 ];

dsh = zeros(npe, dim, npe);
for gp = 1 : npe
   dsh(1,:,gp) = [-1 * (1 - xg(gp,2)) /4 * 2/dx , -1 * (1 - xg(gp,1)) /4 * 2/dy];
   dsh(2,:,gp) = [+1 * (1 - xg(gp,2)) /4 * 2/dx , -1 * (1 + xg(gp,1)) /4 * 2/dy];
   dsh(3,:,gp) = [+1 * (1 + xg(gp,2)) /4 * 2/dx , +1 * (1 + xg(gp,1)) /4 * 2/dy];
   dsh(4,:,gp) = [-1 * (1 + xg(gp,2)) /4 * 2/dx , +1 * (1 - xg(gp,1)) /4 * 2/dy];
end

b_mat = zeros(nvoi, npe*dim, npe);
for gp = 1 : npe
  for i = 1 : npe
      b_mat(1, [i*dim - 1, i*dim + 0], gp) = [dsh(i, 1, gp), 0            ];
      b_mat(2, [i*dim - 1, i*dim + 0], gp) = [0            , dsh(i, 2, gp)];
      b_mat(3, [i*dim - 1, i*dim + 0], gp) = [dsh(i, 2, gp), dsh(i, 1, gp)];
  end
end

endfunction

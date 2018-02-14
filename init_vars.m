function init_vars()

global elements
global coordinates
global bc_nods
global bc_y0
global bc_y1
global bc_x0
global bc_x1
global bc_y0_per
global bc_y1_per
global xg
global wg
global b_mat
global dsh
global npe
global dim
global nvoi

global nx
global ny
global lx
global ly
global dx
global dy

global X0Y0_nod
global X1Y0_nod
global X1Y1_nod
global X0Y1_nod

elements = zeros((nx-1)*(ny-1), npe);
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

bc_y0 = [1 : 1 : nx];
bc_y1 = [(ny-1)*nx + 1 : 1 : nx*ny];
bc_y0_per = [2 : 1 : nx-1];
bc_y1_per = [(ny-1)*nx + 2 : 1 : nx*ny-1];
bc_x0 = [nx + 1 : nx : (ny-2)*nx + 1];
bc_x1 = [2*nx : nx : (ny-1)*nx];
bc_nods = [bc_y0, bc_y1, bc_x0, bc_x1]';

coordinates = zeros(nx*ny, dim);
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

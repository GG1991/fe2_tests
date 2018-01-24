function init_vars()

global elements
global coordinates
global bc_nods
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

elements = zeros((nx-1)*(ny-1), npe);
for i = 1 : (ny-1)
  for j = 1 : (nx-1)
    elements((i-1)*(nx-1) + j, 1) = j    + (i-1)*nx;
    elements((i-1)*(nx-1) + j, 2) = j+1  + (i-1)*nx;
    elements((i-1)*(nx-1) + j, 3) = j+1  + i*nx;
    elements((i-1)*(nx-1) + j, 4) = j    + i*nx;
  end
end

coordinates = zeros(nx*ny, dim);
for i = 1 : ny
  for j = 1 : nx
    coordinates((i - 1)*nx + j, 1) = (j - 1)*dx;
    coordinates((i - 1)*nx + j, 2) = (i - 1)*dy;
  end
end

wg = [0.25 , 0.25, 0.25 , 0.25];

xg = [ -0.577350269189626, -0.577350269189626;
       +0.577350269189626, -0.577350269189626;
       +0.577350269189626, +0.577350269189626;
       -0.577350269189626, +0.577350269189626 ];

dsh = zeros(npe, dim, npe);
for gp = 1 : npe
   dsh(1,1,gp) = -1 * (1 - xg(gp,2)) /4 * 2/dx;
   dsh(2,1,gp) = +1 * (1 - xg(gp,2)) /4 * 2/dx;
   dsh(3,1,gp) = +1 * (1 + xg(gp,2)) /4 * 2/dx;
   dsh(4,1,gp) = -1 * (1 + xg(gp,2)) /4 * 2/dx;
   dsh(1,2,gp) = -1 * (1 - xg(gp,1)) /4 * 2/dy;
   dsh(2,2,gp) = -1 * (1 + xg(gp,1)) /4 * 2/dy;
   dsh(3,2,gp) = +1 * (1 + xg(gp,1)) /4 * 2/dy;
   dsh(4,2,gp) = +1 * (1 - xg(gp,1)) /4 * 2/dy;
end

b_mat = zeros(nvoi, npe*dim, npe);
for gp = 1 : npe
  for i = 1 : npe
      b_mat(1, i*dim - 1, gp) = dsh(i, 1, gp);
      b_mat(1, i*dim + 0, gp) = 0;
      b_mat(2, i*dim - 1, gp) = 0;
      b_mat(2, i*dim + 0, gp) = dsh(i, 2, gp);
      b_mat(3, i*dim - 1, gp) = dsh(i, 2, gp);
      b_mat(3, i*dim + 0, gp) = dsh(i, 1, gp);
  end
end

endfunction

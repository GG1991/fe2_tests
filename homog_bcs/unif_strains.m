nx = 10;
ny = 10;
lx = 3;
ly = 3;
dx = lx / (nx - 1);
dy = ly / (ny - 1);
npe = 4;
dim = 2;

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

#elements
#coordinates

u_n = zeros(nx*ny*dim, 1);

bc_y0 = [1 : 1 : nx];
bc_y1 = [(ny-1)*nx + 1 : 1 : nx*ny];
bc_x0 = [nx + 1 : nx : (ny-2)*nx + 1];
bc_x1 = [2*nx : nx : (ny-1)*nx];
bc_index = [bc_y0, bc_y1, bc_x0, bc_x1];

dir_n = zeros(nx*ny*dim, 1);

[jac, res] = ass_unifstrains (elements, coordinates, nx, ny, lx, ly, u_n);

%figure();
%spy(jac)
%print -djpg filename.jpg 

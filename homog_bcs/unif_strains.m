nx = 3;
ny = 3;
dx = 2;
dy = 2;
npe = 4;

elements = zeros((nx-1)*(ny-1), npe);
for i = 1 : (ny-1)
  for j = 1 : (nx-1)
    elements((i-1)*(nx-1) + j, 1) = j    + (i-1)*nx;
    elements((i-1)*(nx-1) + j, 2) = j+1  + (i-1)*nx;
    elements((i-1)*(nx-1) + j, 3) = j+1  + i*nx;
    elements((i-1)*(nx-1) + j, 4) = j    + i*nx;
  end
end

elements

[jac, res] = ass_unifstrains (elements, nx, ny, dx, dy);

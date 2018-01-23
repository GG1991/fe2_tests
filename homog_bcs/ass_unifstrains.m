%  ensambla el Jacobiano (J) y el residuo (r)

function [jac, res] = ass_unifstrains (elements, coordinates, strain_mac, bc_nods, nx, ny, lx, ly, u_n)

nnods = nx*ny;
dim = 2;
npe = size(elements, 2);
nvoi = 3;
dx = lx / (nx - 1);
dy = ly / (ny - 1);

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

jac = sparse(dim*nnods, dim*nnods);
res = zeros(dim*nnods, 1);

for e = 1 : size(elements, 1) 

    u_e = u_n([elements(e, :)*dim - 1, elements(e, :)*dim]);
    [ejac, eres] = elemental (e, elements, coordinates, lx, ly, xg, wg, dsh, b_mat, u_e);

    jac([elements(e,:)*dim - 1, elements(e,:)*dim], [elements(e,:)*dim - 1, elements(e,:)*dim]) += ejac;
    res([elements(e,:)*dim - 1, elements(e,:)*dim]) += eres;

end

u_d = zeros(size(bc_nods, 1)*dim);
for n = 1 : size(bc_nods, 1) 
end

endfunction

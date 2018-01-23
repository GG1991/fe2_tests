%  ensambla el Jacobiano (J) y el residuo (r)

function [jac, res] = ass_unifstrains (elements, nx, ny, dx, dy)

nnods = nx*ny;
dim = 2;
npe = size(elements, 2);

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

jac = sparse(dim*nnods, dim*nnods);
res = zeros(dim*nnods, 1);

for e = 1 : size(elements, 1) 

    [ejac, eres] = calc_elemental (e, xg, wg, dsh);

    jac([elements(e,:)*dim,(elements(e,:)*dim) - 1], [elements(e,:)*dim,(elements(e,:)*dim)]) += ejac;
    res([elements(e,:)*dim,(elements(e,:)*dim) - 1]) += eres;

end

endfunction

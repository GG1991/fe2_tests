% Steffiness matrix
function [jac, res] = ass_steff (u)

global dim
global npe
global nn
global nelem
global elements

jac = sparse(nn*dim, nn*dim);
res = zeros(nn*dim, 1);
u_e = zeros(npe*dim, 1);
ind = zeros(npe*dim, 1);

for e = 1 : nelem 
    u_e([1:2:npe*dim, 2:2:npe*dim]) = u([elements(e, :)*dim - 1, elements(e, :)*dim + 0]);

    [jac_e, res_e] = elemental (e, u_e);
    ind = [elements(e,:)*dim - 1; elements(e,:)*dim - 0](:);

    jac(ind, ind) += jac_e;
    res(ind) += res_e;
end

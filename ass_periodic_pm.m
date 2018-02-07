% Penalty method
function [K, f] = ass_periodic_pm (strain_mac, u_n)

global elements
global coordinates
global stress
global strain
global bc_nods
global bc_y0_per
global bc_y1_per
global bc_x0
global bc_x1
global lx
global ly
global dx
global dy
global nx
global ny
global nnods
global npe
global dim
global nelem
global size_tot

w = 1.0e10;

K = sparse(size_tot, size_tot);
f = zeros(size_tot, 1);
u_e = zeros(npe*dim, 1);
ind = zeros(npe*dim, 1);

for e = 1 : nelem 

    u_e([1:2:npe*dim]) = u_n([elements(e, :)*dim - 1]); %set x vals
    u_e([2:2:npe*dim]) = u_n([elements(e, :)*dim + 0]); %set y vals

    [K_e, f_e] = elemental (e, u_e);
    ind = [elements(e,:)*dim - 1; elements(e,:)*dim - 0](:);

    K(ind, ind) += K_e;

end

X0Y0_nod = 1;
X1Y0_nod = nx;
X1Y1_nod = nx*ny;
X0Y1_nod = (ny-1)*nx + 1;

u_dif_y0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [0.0, ly]';
u_dif_x0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [lx , 0.0]';

u_X0Y0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [0.0, 0.0]';
u_X1Y0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [lx , 0.0]';
u_X1Y1 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [lx , ly ]';
u_X0Y1 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [0.0, ly ]';

% fb +- µ
f(bc_y0_per*dim - 1) -= w*u_dif_y0(1); %x
f(bc_y1_per*dim - 1) += w*u_dif_y0(1); %x
f(bc_y0_per*dim - 0) -= w*u_dif_y0(2); %y
f(bc_y1_per*dim - 0) += w*u_dif_y0(2); %y

f(bc_x0*dim - 1) -= w*u_dif_x0(1); %x
f(bc_x1*dim - 1) += w*u_dif_x0(1); %x
f(bc_x0*dim - 0) -= w*u_dif_x0(2); %y
f(bc_x1*dim - 0) += w*u_dif_x0(2); %y

f([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) = u_X0Y0; % x & y
f([X1Y0_nod*dim - 1, X1Y0_nod*dim + 0]) = u_X1Y0; % x & y
f([X1Y1_nod*dim - 1, X1Y1_nod*dim + 0]) = u_X1Y1; % x & y
f([X0Y1_nod*dim - 1, X0Y1_nod*dim + 0]) = u_X0Y1; % x & y

for n = 1 : max(size(bc_y0_per))
  for d = 0 : 1
    K(bc_y1_per(n)*dim - d, bc_y1_per(n)*dim - d) += w; %diag
    K(bc_y1_per(n)*dim - d, bc_y0_per(n)*dim - d) -= w;
    K(bc_y0_per(n)*dim - d, bc_y0_per(n)*dim - d) += w; %diag
    K(bc_y0_per(n)*dim - d, bc_y1_per(n)*dim - d) -= w;
  end
end

for n = 1 : max(size(bc_x0))
  for d = 0 : 1
    K(bc_x1(n)*dim - d, bc_x1(n)*dim - d) += w; %diag
    K(bc_x1(n)*dim - d, bc_x0(n)*dim - d) -= w;
    K(bc_x0(n)*dim - d, bc_x0(n)*dim - d) += w; %diag
    K(bc_x0(n)*dim - d, bc_x1(n)*dim - d) -= w;
  end
end

K([X0Y0_nod*dim - 1; X0Y0_nod*dim - 0], :) = 0.0;
K([X1Y0_nod*dim - 1; X1Y0_nod*dim - 0], :) = 0.0;
K([X1Y1_nod*dim - 1; X1Y1_nod*dim - 0], :) = 0.0;
K([X0Y1_nod*dim - 1; X0Y1_nod*dim - 0], :) = 0.0;
for d = 0 : 1
 K(X0Y0_nod*dim - d, X0Y0_nod*dim - d) = 1.0;
 K(X1Y0_nod*dim - d, X1Y0_nod*dim - d) = 1.0;
 K(X1Y1_nod*dim - d, X1Y1_nod*dim - d) = 1.0;
 K(X0Y1_nod*dim - d, X0Y1_nod*dim - d) = 1.0;
end

endfunction

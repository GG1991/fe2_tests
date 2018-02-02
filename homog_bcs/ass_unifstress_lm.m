% Lagrange multiplier method
function [jac, res] = ass_unifstress_lm (strain_mac, u_n)

global elements
global coordinates
global stress
global strain
global bc_nods
global bc_y0_per
global bc_y1_per
global bc_y0
global bc_y1
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

global lam_1
global lam_2
global lam_3

area = lx*ly;

jac = sparse(size_tot, size_tot);
res = zeros(size_tot, 1);
u_e = zeros(npe*dim, 1);
ind = zeros(npe*dim, 1);

for e = 1 : nelem 

    u_e([1:2:npe*dim]) = u_n([elements(e, :)*dim - 1]); %set x vals
    u_e([2:2:npe*dim]) = u_n([elements(e, :)*dim + 0]); %set y vals

    [jac_e, res_e] = elemental (e, u_e);
    for n = 1 : npe 
      for d = 0 : 1
        ind(n*dim - d) = elements(e,n)*dim - d;
      end
    end

    jac(ind, ind) += jac_e;
    res(ind) += res_e;

end

X0Y0_nod = 1;
X1Y0_nod = nx;
X1Y1_nod = nx*ny;
X0Y1_nod = (ny-1)*nx + 1;

u_X0Y0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [0.0, 0.0]';
u_X1Y0 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [lx , 0.0]';
u_X1Y1 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [lx , ly ]';
u_X0Y1 = [strain_mac(1) strain_mac(3)/2 ; strain_mac(3)/2 strain_mac(2)] * [0.0, ly ]';

ax = dx/area;
ay = dy/area;

lam_1 = u_n(nx*ny*dim + 1);
lam_2 = u_n(nx*ny*dim + 2);
lam_3 = u_n(nx*ny*dim + 3);

res(nx*ny*dim + 1)  = + sum(u_n(bc_x1*dim     - 1))*ay   - sum(u_n(bc_x0*dim     - 1))*ay   - strain_mac(1); % ux x nx - exx
res(nx*ny*dim + 2)  = + sum(u_n(bc_y1_per*dim - 0))*ax   - sum(u_n(bc_y0_per*dim - 0))*ax   - strain_mac(2); % uy x ny - eyy
res(nx*ny*dim + 3)  = + sum(u_n(bc_y1_per*dim - 1))*ax/2 - sum(u_n(bc_y0_per*dim - 1))*ax/2;                 % 1/2 ux x ny
res(nx*ny*dim + 3) += + sum(u_n(bc_x1*dim     - 0))*ay/2 - sum(u_n(bc_x0*dim     - 0))*ay/2;                 % 1/2 uy x nx
res(nx*ny*dim + 3) += - strain_mac(3)/2; % - exy

% corners
res(nx*ny*dim + 1) += + u_n(X1Y0_nod*dim - 1)*ay/2 + u_n(X1Y1_nod*dim - 1)*ay/2 - u_n(X0Y1_nod*dim -1)*ay/2; % X0 X1 faces
res(nx*ny*dim + 2) += - u_n(X1Y0_nod*dim - 0)*ax/2 + u_n(X1Y1_nod*dim - 0)*ax/2 + u_n(X0Y1_nod*dim -0)*ax/2; % Y0 Y1 faces
res(nx*ny*dim + 3) += - u_n(X1Y0_nod*dim - 1)*ax/4 + u_n(X1Y1_nod*dim - 1)*ax/4 + u_n(X0Y1_nod*dim -1)*ax/4; % Y0 Y1 faces
res(nx*ny*dim + 3) += + u_n(X1Y0_nod*dim - 0)*ay/4 + u_n(X1Y1_nod*dim - 0)*ay/4 - u_n(X0Y1_nod*dim -0)*ay/4; % X0 X1 faces

res(bc_y1_per*dim - 1) -= +ax/2*lam_3;
res(bc_y1_per*dim - 0) -= +ax  *lam_2;
res(bc_y0_per*dim - 1) -= -ax/2*lam_3;
res(bc_y0_per*dim - 0) -= -ax  *lam_2;

res(bc_x1*dim     - 1) -= +ay  *lam_1;
res(bc_x1*dim     - 0) -= +ay/2*lam_3;
res(bc_x0*dim     - 1) -= -ay  *lam_1;
res(bc_x0*dim     - 0) -= -ay/2*lam_3;

% corners
res([X1Y0_nod X1Y1_nod X0Y1_nod]*dim - 1) += -[+ay/2 +ay/2 -ay/2]'*lam_1; % X0 X1 faces
res([X1Y0_nod X1Y1_nod X0Y1_nod]*dim - 0) += -[-ax/2 +ax/2 +ax/2]'*lam_2; % Y0 Y1 faces
res([X1Y0_nod X1Y1_nod X0Y1_nod]*dim - 1) += -[-ax/4 +ax/4 +ax/4]'*lam_3; % Y0 Y1 faces
res([X1Y0_nod X1Y1_nod X0Y1_nod]*dim - 0) += -[+ay/4 +ay/4 -ay/4]'*lam_3; % X0 X1 faces

res([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) = u_n([X0Y0_nod*dim - 1, X0Y0_nod*dim + 0]) - u_X0Y0; % x & y

jac(nx*ny*dim + 1, bc_x1*dim     - 1) = +ay;
jac(nx*ny*dim + 1, bc_x0*dim     - 1) = -ay;
jac(nx*ny*dim + 2, bc_y1_per*dim - 0) = +ax;
jac(nx*ny*dim + 2, bc_y0_per*dim - 0) = -ax;
jac(nx*ny*dim + 3, bc_x1*dim     - 0) = +ay/2;
jac(nx*ny*dim + 3, bc_x0*dim     - 0) = -ay/2;
jac(nx*ny*dim + 3, bc_y1_per*dim - 1) = +ax/2;
jac(nx*ny*dim + 3, bc_y0_per*dim - 1) = -ax/2;

% corners
jac(nx*ny*dim + 1,[X1Y0_nod X1Y1_nod X0Y1_nod]*dim - 1) += [+ay/2 +ay/2 -ay/2]; % X0 X1 faces
jac(nx*ny*dim + 2,[X1Y0_nod X1Y1_nod X0Y1_nod]*dim - 0) += [-ax/2 +ax/2 +ax/2]; % Y0 Y1 faces
jac(nx*ny*dim + 3,[X1Y0_nod X1Y1_nod X0Y1_nod]*dim - 1) += [-ax/4 +ax/4 +ax/4]; % Y0 Y1 faces
jac(nx*ny*dim + 3,[X1Y0_nod X1Y1_nod X0Y1_nod]*dim - 0) += [+ay/4 +ay/4 -ay/4]; % X0 X1 faces

for n = 1 : max(size(bc_y0_per))
  jac(bc_y1_per(n)*dim - 1, nx*ny*dim + [1 2 3]) = -[0   0    +ax/2];
  jac(bc_y1_per(n)*dim - 0, nx*ny*dim + [1 2 3]) = -[0   +ax  0    ];
  jac(bc_y0_per(n)*dim - 1, nx*ny*dim + [1 2 3]) = -[0   0    -ax/2];
  jac(bc_y0_per(n)*dim - 0, nx*ny*dim + [1 2 3]) = -[0   -ax  0    ];

  jac(bc_x1(n)*dim - 1    , nx*ny*dim + [1 2 3]) = -[+ay 0    0    ];
  jac(bc_x1(n)*dim - 0    , nx*ny*dim + [1 2 3]) = -[0   0    +ay/2];
  jac(bc_x0(n)*dim - 1    , nx*ny*dim + [1 2 3]) = -[-ay 0    0    ];
  jac(bc_x0(n)*dim - 0    , nx*ny*dim + [1 2 3]) = -[0   0    -ay/2];
end

% corners
jac(X1Y0_nod*dim - 1, nx*ny*dim + [1 2 3]) += -[+ay/2   0.0   -ax/4]; % X1 face x dir
jac(X1Y0_nod*dim - 0, nx*ny*dim + [1 2 3]) += -[ 0.0   -ax/2  +ay/4]; % X1 face y dir
jac(X1Y1_nod*dim - 1, nx*ny*dim + [1 2 3]) += -[+ay/2   0.0   +ax/4]; % X1 face x dir
jac(X1Y1_nod*dim - 0, nx*ny*dim + [1 2 3]) += -[ 0.0   +ax/2  +ay/4]; % X1 face y dir

jac(X0Y1_nod*dim - 1, nx*ny*dim + [1 2 3]) += -[-ay/2   0.0   +ax/4]; % X0 face x dir
jac(X0Y1_nod*dim - 0, nx*ny*dim + [1 2 3]) += -[ 0.0   +ax/2  -ay/4]; % X0 face y dir

jac([X0Y0_nod*dim - 1; X0Y0_nod*dim - 0], :) = 0.0;
for d = 0 : 1
 jac(X0Y0_nod*dim - d, X0Y0_nod*dim - d) = 1.0;
end

endfunction

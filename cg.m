function x_1 = cg(A, b, x_0)

K_inv = diag(1./diag(A));

tol_min = 1.0e-5;
its_max = 10;

its = 1;
tol = 1.0;

r_0 = A*x_0 - b;
while (its < its_max && tol > tol_min)

  rho_1 = r_0'*K_inv*r_0;
  if (its == 1)
    p_1 = K_inv*r_0;
  else
    p_1 = K_inv*r_0 + (rho_1/rho_0)*p_0;
  endif

  d_1 = rho_1 / (p_1'*A*p_1);

  x_1 = x_0 - d_1*p_1;
  r_1 = r_0 - d_1*A*p_1;
  x_0 = x_1;
  p_0 = p_1;
  r_0 = r_1;
  rho_0 = rho_1;

  its += 1;
  tol = norm(A*x_1 - b);
  printf("%e\n",tol);

endwhile
printf("%d\n",its);

endfunction

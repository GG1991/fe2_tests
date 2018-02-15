function [x_1, err, its] = cg(A, b, x_0, min_tol, max_its)

K_inv = diag(1./diag(A));

its = 1;
err = 1.0;

r_0 = A*x_0 - b;
while (its < max_its && err > min_tol)

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
  err = norm(A*x_1 - b);

endwhile

endfunction

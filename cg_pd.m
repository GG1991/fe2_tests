function [x, err, its] = cg_pd(A, b, x, min_tol, max_its)

its = 1;
M_i = diag(1./diag(A));

r_0 = b - A*x;
z_0 = M_i*r_0; 
err = norm(r_0);
p = z_0;
while (its < max_its && err > min_tol)
  Ap = A*p;
  alpha = (r_0' * z_0) / (p' * Ap);
  x = x + alpha*p;
  r_1 = r_0 - alpha*Ap;
  z_1 = M_i*r_1; 
  beta = (z_1' * r_1)/(z_0' * r_0);
  p = z_1 + beta * p;
  r_0 = r_1;
  z_0 = z_1;
  its += 1;
  err = norm(r_1);
endwhile

endfunction

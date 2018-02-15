function [x, err, its] = cg(A, b, x, min_tol, max_its)

its = 1;

r_0 = b - A*x;
err = norm(r_0);
p = r_0;
while (its < max_its && err > min_tol)
  Ap = A*p;
  alpha = (r_0' * r_0) / (p' * Ap);
  x = x + alpha*p;
  r_1 = r_0 - alpha*Ap;
  beta = (r_1' * r_1)/(r_0' * r_0);
  p = r_1 + beta * p;
  r_0 = r_1;
  its += 1;
  err = norm(r_1);
endwhile

endfunction

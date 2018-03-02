function [x1, x2, err, its] = cg_usawa(A, B, b1, b2, x1, x2, min_tol, max_its)

% solves
%
% |A B||x1| = |b1|
% |B 0||x2|   |b2|
%
% https://en.wikipedia.org/wiki/Uzawa_iteration

its = 1;

r2 = B'*x1 - b2;
err = norm(r2);
p2 = r2;
while (its < max_its && err > min_tol)

  p1 = A\(B*p2);
  a2 = B'*p1;
  alpha = (p2' * r2) / (p2' * a2);
  x2 = x2 + alpha*p2;
  x1 = x1 - alpha*p1;
  r2 = r2 - alpha*a2;

  beta = (r2' * a2)/(p2' * a2);
  p2 = r2 - beta * p2;

  its += 1;
  err = norm(r2);
endwhile

endfunction

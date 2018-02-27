function [u, tol, its] = cg_pgs(A, f, u_s, min_tol, max_its)

D = diag(diag(A));
C1 = tril(A);
C2 = inv(D)*triu(A);
	
u = u_s;
r = f - A * u;
p = C2 \ (C1 \ r);
norm_f = norm(f);

tol = norm(r);
its = 0;
while((tol > min_tol) && (its < max_its))
  a = A * p;
  a_dot_p = a' * p;
  lambda = (r' * p) / a_dot_p;
  u = u + lambda * p;
  r = r - lambda * a;
  inv_C_times_r = C2 \ (C1 \ r);
  p = inv_C_times_r - ((inv_C_times_r' * a) / a_dot_p) * p;
  tol = norm(r);
  its += 1;
end

function [x, tol, its] = jacobi(A,b,x0)

D = diag(diag(A));
N = D - A;

tol_min = 1.0e-5;
its_max = 1000;

its = 1;
tol = 100;

x_old = x0;
while (its < its_max && tol > tol_min)
  x = D \ (N*x_old + b);
  x_old = x; 
  tol = norm(A*x - b);
  its += 1;
end

endfunction

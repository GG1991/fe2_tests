n = 40;
A = sparse(n,n);
A_e = [1 -1; -1 1]; 

for i = 1 : n - 1
A([i i+1], [i i+1]) += A_e;
end
A(1,:) = 0;
A(:,1) = 0;
A(1,1) = 1;
A(n,:) = 0;
A(:,n) = 0;
A(n,n) = 1;
%full(A)
b = ones(n,1);
b(1) = 0;
b(n) = 0;

%x = A \ b

x = zeros(n,1);
tic();
[x, tol, its] = jacobi(A,b,x);
time = toc();
printf("\033[31m jacobi\n", tol);
printf("\033[32m tol :%e\n", tol);
printf("\033[32m its :%d\n", its);
printf("\033[32m time:%f\n", time);
save -ascii sol_jacobi.dat x
printf("\n");

x = zeros(n,1);
tic();
[x, tol, its] = cg(A,b,x);
time = toc();
printf("\033[31m cg\n", tol);
printf("\033[32m tol :%e\n", tol);
printf("\033[32m its :%d\n", its);
printf("\033[32m time:%f\n", time);
save -ascii sol_cg.dat x

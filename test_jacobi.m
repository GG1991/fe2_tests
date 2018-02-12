n = 6;
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
full(A)
b = ones(n,1);
b(1) = 0;
b(n) = 0;

x = A \ b

x1 = zeros(n,1);
x1 = jacobi(A,b,x1)

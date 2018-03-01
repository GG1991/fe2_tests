function [x] = my_lu(A, b)

[L, U, p] = lu(A);

%foward
y = zeros(size(A,1),1);
for i = 1 : size(A,1) 
  y(i) = b(i);
  for j = 1 : (i-1)
    y(i) -= L(i,j)*y(j);
  end
end

%backward
x = zeros(size(A,1),1);
for i = size(A,1) : -1 : 1
  x(i) = y(i);
  for j = i : (size(A,1)-1)
    x(i) -= U(i,j)*x(j);
  end
  x(i) = x(i) / U(i,i);
end


endfunction

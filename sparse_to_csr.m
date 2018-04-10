function [p, i, x] = sparse_to_csr(A)
  [i, j, x] = find(A);
  p = cumsum([1 ; accumarray(j, 1)]);
endfunction 

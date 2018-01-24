function [strain_ave, stress_ave] = average ()

global nelem
global strain
global stress

strain_ave = [0 0 0];
stress_ave = [0 0 0];

for e = 1 : nelem

  strain_ave += strain(e,:);
  stress_ave += stress(e,:);

end

strain_ave /= nelem;
stress_ave /= nelem;

endfunction

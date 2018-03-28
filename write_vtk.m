function write_vtk(file_name, u_vec)

global elements; global coordinates; global elem_type; global nn; global nelem
global strain; global stress; global res; global dim; global int_vars; global wg

fm = fopen (file_name, "w");
fprintf(fm, '# vtk DataFile Version 2.0\n');
fprintf(fm, 'fe2_tests\n');
fprintf(fm, 'ASCII\n');

fprintf(fm,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fm,'POINTS %d FLOAT\n', nn);
for n = 1 : nn
   fprintf(fm,'%f %f 0.0\n',coordinates(n,:));
end
fprintf(fm,'CELLS %d %d\n', nelem, 5*nelem);
for e = 1 : nelem
   fprintf(fm,'4 %d %d %d %d\n',elements(e,:)-1);
end
fprintf(fm,'CELL_TYPES %d\n', nelem);
for e = 1 : nelem
   fprintf(fm,'9\n');
end

fprintf(fm, 'POINT_DATA %d\n', nn);
fprintf(fm, 'VECTORS disp FLOAT\n');
%fprintf(fm, 'LOOKUP_TABLE default\n');
for n = 1 : nn
   fprintf(fm,'%f %f 0.0\n',u_vec([n*dim - 1, n*dim + 0]));
end
%fprintf(fm, 'VECTORS residue FLOAT\n');
%%fprintf(fm, 'LOOKUP_TABLE default\n');
%for n = 1 : nn
%   fprintf(fm,'%f %f 0.0\n',res([n*dim - 1, n*dim + 0]));
%end

fprintf(fm, 'CELL_DATA %d\n', nelem);
fprintf(fm, 'SCALARS type FLOAT\n');
fprintf(fm, 'LOOKUP_TABLE default\n');
for e = 1 : nelem
   fprintf(fm,'%f\n',elem_type(e));
end

fprintf(fm, 'SCALARS hard FLOAT\n');
fprintf(fm, 'LOOKUP_TABLE default\n');
for e = 1 : nelem
  eps_p_ave = [0 0 0];
  alpha_ave = 0;
  for gp = 1 : 4
    eps_p_ave += int_vars((e-1)*4 + gp,[1 2 3]) * wg(gp);
    alpha_ave += int_vars((e-1)*4 + gp,[7]) * wg(gp);
  endfor
  eps_p_ave /= sum(wg);
  alpha_ave /= sum(wg);
  fprintf(fm,'%f\n', norm(eps_p_ave));
endfor

fprintf(fm, 'TENSORS strain FLOAT\n');
for e = 1 : nelem
   fprintf(fm,'%f %f 0 %f %f 0 0 0 0\n',strain(e,1),strain(e,3),strain(e,3),strain(e,2));
end

fprintf(fm, 'TENSORS stress FLOAT\n');
for e = 1 : nelem
   fprintf(fm,'%f %f 0 %f %f 0 0 0 0\n',stress(e,1),stress(e,3),stress(e,3),stress(e,2));
end

fclose(fm);

endfunction

function write_vtk(file_name, u_vec)

global elements
global coordinates
global elem_type
global nnods
global nelem
global strain
global stress
global dim

fm = fopen (file_name, "w");
fprintf(fm, '# vtk DataFile Version 2.0\n');
fprintf(fm, 'fe2_tests\n');
fprintf(fm, 'ASCII\n');

fprintf(fm,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fm,'POINTS %d FLOAT\n', nnods);
for n = 1 : nnods
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

fprintf(fm, 'POINT_DATA %d\n', nnods);
fprintf(fm, 'VECTORS disp FLOAT\n');
%fprintf(fm, 'LOOKUP_TABLE default\n');
for n = 1 : nnods
   fprintf(fm,'%f %f 0.0\n',u_vec([n*dim - 1, n*dim + 0]));
end

fprintf(fm, 'CELL_DATA %d\n', nelem);
fprintf(fm, 'SCALARS type FLOAT\n');
fprintf(fm, 'LOOKUP_TABLE default\n');
for e = 1 : nelem
   fprintf(fm,'%f\n',elem_type(e));
end

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

function dist_val = distance (e, elements, coordinates, lx, ly)

global elements
global coordinates
global lx
global ly

center = [lx / 2-0.0, ly / 2+0.0];
centroid = [0.0, 0.0];

for i = 1 : size(elements, 2)
  centroid += coordinates(elements(e,i),:);
end
centroid /= size(elements, 2);

dist_val = norm(centroid - center);

endfunction

function [] = read_comm_line(nargin, argv)

global nx; global ny; global dim; global nvoi; global bc_y0; global bc_x0;
global solver; global bc_type; global nexp; global size_tot;

for i = 1 : nargin
 if (strcmp(argv(){i}, "-cg"))
  solver = "cg";
 elseif (strcmp(argv(){i}, "-cg_pd"))
  solver = "cg_pd";
 elseif (strcmp(argv(){i}, "-cg_pgs"))
  solver = "cg_pgs";
 elseif (strcmp(argv(){i}, "-lu"))
  solver = "lu";
 elseif (strcmp(argv(){i}, "-my_lu"))
  solver = "my_lu";
 elseif (strcmp(argv(){i}, "-ustrain"))
  bc_type = "ustrain";
 elseif (strcmp(argv(){i}, "-ustress"))
  bc_type = "ustress";
 elseif (strcmp(argv(){i}, "-per_lm"))
  bc_type = "per_lm";
 elseif (strcmp(argv(){i}, "-per_ms"))
  bc_type = "per_ms";
 elseif (strcmp(argv(){i}, "-nx"))
  nx = str2num(argv(){i+1});
 elseif (strcmp(argv(){i}, "-ny"))
  ny = str2num(argv(){i+1});
 elseif (strcmp(argv(){i}, "-nexp"))
  nexp = str2num(argv(){i+1});
 endif
end

endfunction

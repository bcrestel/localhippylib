fid = fopen('arolla.geo','w');

TT = load('arolla100_modified.dat');
points = [TT(:,1:2); TT(end-1:-1:2,[1,3])];
npoints = size(points,1);

id = [1:npoints]';

fprintf(fid, 'lc = %.1f;\n', 25); 

tmp = [ id, points, zeros(npoints,1)];
fprintf(fid, 'Point(%d) = {%.1f, %.1f, %.1f, lc};\n', tmp');
nseg = npoints;

tmp = [id, id, id([2:end, 1])];
fprintf(fid, 'Line(%d) = {%d, %d};\n', tmp');
fprintf(fid, 'Line Loop(%d) = {',nseg+1);
fprintf(fid, '%d, ', id(1:end-1));
fprintf(fid, '%d};\n', id(end));
fprintf(fid, 'Plane Surface(%d) = {%d};\n', nseg+2, nseg+1);

fprintf(fid, '\nPhysical Line(1) = {');
fprintf(fid, '%d, ', id(1:end/2-1));
fprintf(fid, '%d};\n', id(end/2));

fprintf(fid, '\nPhysical Line(2) = {');
fprintf(fid, '%d, ', id(end/2+1:end-1));
fprintf(fid, '%d};\n', id(end));

fprintf(fid, '\n Physical Surface(3) = {%d};', nseg+2);





fclose(fid);
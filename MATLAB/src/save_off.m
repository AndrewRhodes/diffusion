


function save_off(vertices, faces, FILE_NAME)

NumVertices = length(vertices);
NumFaces = length(faces);


h = waitbar(1, sprintf('Saving as .off file...'));
pause(1);

% pathName = '3D_Model';


fid = fopen(fullfile(FILE_NAME), 'w+');
fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d %d\n', NumVertices, NumFaces, 0);
fprintf(fid, '%6f %6f %6f\n', vertices');
fprintf(fid, '3 %d %d %d\n', faces'-1);
fclose(fid);

close(h);


end




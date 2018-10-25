


function write_pointcloud_ply(filename, verts)




fileID = fopen(filename,'w');


fprintf(fileID, ...
    ['ply\n', ...
    'format ascii 1.0\n', ...
    'element vertex %u\n', ...
    'property float32 x\n', ...
    'property float32 y\n', ...
    'property float32 z\n', ...
    'end_header\n'], ...
    size(verts,1));


for i= 1 : size(verts,1)
    fprintf(fileID, ...
        ['%.6f ', ...
        '%.6f ', ...
        '%.6f\n'], ...
        verts(i,1),verts(i,2),verts(i,3));
end


fclose(fileID);


end







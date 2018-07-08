

close all
clear 
clc

global ProjectRoot; % Additional Paths


addpath('~/GitProjects/matlab-utilities/')
addpath('~/GitProjects/mesh_resampling_toolbox/MATLAB_Modules/')
% addpath(genpath('~/Documents/Software/cp_matrices/'))
addpath('~/Desktop/MLIDAR-1.0/MATLAB_Modules/')
addpath(genpath('~/GitProjects/pose/MLIDAR-2.2/'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear PointCloud
DirectoryName = 'dragon';
ModelName = 'Dragon_e1';
DownSampleFacesNum = 50000;

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(DirectoryName,'/',ModelName,'.ply');

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
PointCloud = findMeshResolution(PointCloud, 'Model');

% % % PointCloud.Location = PointCloud.Location .* (1/PointCloud.Resolution);

PC.faces = PointCloud.Face;
PC.vertices = PointCloud.Location;

PC2 = reducepatch(PC, DownSampleFacesNum);

clear PointCloud

PointCloud.Face = PC2.faces;
PointCloud.Location = PC2.vertices;

[VerticesOut, FacesOut] = clearMeshDuplicates(PointCloud.Location, PointCloud.Face );

clear PointCloud

PointCloud.Face = FacesOut;
PointCloud.Location = VerticesOut;
PointCloud = findMeshNormals(PointCloud);
PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud = findMeshResolution(PointCloud, 'Model');
PointCloud


ply_write(strcat(ModelName,'_',num2str(DownSampleFacesNum),'.ply'), PointCloud.Face, PointCloud.Location)
save_off(PointCloud.Location, PointCloud.Face, strcat(ModelName,'_',num2str(DownSampleFacesNum),'.off'))







% sunvector = [-.5;-.5;1];
sunvector = [0;0;1];
sunvector = sunvector / norm(sunvector);
intensity = 0.9*PointCloud.Normal*sunvector;
intensity(intensity<0) = 0;



figure
colormap('gray')
% trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), 'EdgeColor', 'none');
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), intensity, 'EdgeColor', 'none');
view(0, 90)
axis equal
axis off
hold on













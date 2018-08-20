

global ProjectRoot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear PointCloud
Model = 'bunny/Bunny_e1';

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(Model,'.ply');

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud = findMeshNormals(PointCloud);


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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear PointCloud

Model = 'armadillo/Armadillo_e1_100000';

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(Model,'.ply');

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud = findMeshNormals(PointCloud);


sunvector = [0;0;-1];
sunvector = sunvector / norm(sunvector);
intensity = 0.9*PointCloud.Normal*sunvector;
intensity(intensity<0) = 0;

figure
colormap('gray')
% trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), 'EdgeColor', 'none');
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), intensity, 'EdgeColor', 'none');
view(180, -90)
axis equal
axis off
hold on






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear PointCloud

Model = 'buddha/Buddha_e1';

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(Model,'.ply');

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud = findMeshNormals(PointCloud);
PointCloud = findMeshResolution(PointCloud, 'Model');


sunvector = [0;0;1];
sunvector = sunvector / norm(sunvector);
intensity = 0.9*PointCloud.Normal*sunvector;
intensity(intensity<0) = 0;
intensity(isnan(intensity)) = 0;


figure
colormap('gray')
% trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), 'EdgeColor', 'none');
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), intensity, 'EdgeColor', 'none');
view(-360, 90)
axis equal
axis off
hold on






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear PointCloud

Model = 'dragon/Dragon_e1';

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(Model,'.ply');

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud = findMeshNormals(PointCloud);
% PointCloud = findMeshResolution(PointCloud, 'Model');


sunvector = [0;0;1];
sunvector = sunvector / norm(sunvector);
intensity = 0.9*PointCloud.Normal*sunvector;
intensity(intensity<0) = 0;
intensity(isnan(intensity)) = 0;


figure
colormap('gray')
% trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), 'EdgeColor', 'none');
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), intensity, 'EdgeColor', 'none');
view(-360, 90)
axis equal
axis off
hold on




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear PointCloud

Model = 'itokawa/Itokawa_e1_80000';

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(Model,'.ply');

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud = findMeshNormals(PointCloud);
% PointCloud = findMeshResolution(PointCloud, 'Model');


sunvector = [0;0;1];
sunvector = sunvector / norm(sunvector);
intensity = 0.9*PointCloud.Normal*sunvector;
intensity(intensity<0) = 0;
intensity(isnan(intensity)) = 0;


figure
colormap('gray')
% trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), 'EdgeColor', 'none');
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), intensity, 'EdgeColor', 'none');
view(-200, 90)
axis equal
axis off
hold on
















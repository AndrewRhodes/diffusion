
close all
clear 
clc
global ProjectRoot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear PointCloud
ModelFolder = 'bunny/';
Model = 'Bunny_e1';


FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud = findMeshNormals(PointCloud);

FileLocationNeighbors = strcat('/media/andrew/WDRhodes/diffusiondata/',ModelFolder,'neighbors/');
load(strcat(FileLocationNeighbors, Model,'_Neighbors.mat'), 'Neighbors')

PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);



% sunvector = [-.5;-.5;1];
sunvector = [0;0;1];
sunvector = sunvector / norm(sunvector);
intensity = 0.9*PointCloud.Normal*sunvector;
intensity(intensity<0) = 0;


FileLocation = strcat('/media/andrew/WDRhodes/diffusiondata/',ModelFolder,'SignalNoise/Umbrella_te8/');
% FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'VertexNoise/Mesh_t2_geo/');
load(strcat(FileLocation,'NMSKeypoint.mat'),'NMSKeypoint')



figure
colormap('gray')
% trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), 'EdgeColor', 'none');
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), intensity, 'EdgeColor', 'none');
view(0, 90)
axis equal
axis off
hold on



ScaleLogic = NMSKeypoint.Scale < ( 3 * PointCloud.ResolutionLocal(NMSKeypoint.LocationIndex) );

NMSKeypoint = rmfield(NMSKeypoint, 'Count');
FNames = fieldnames(NMSKeypoint);
for jj = 1 : length(FNames)
    NMSKeypoint.(FNames{jj})(ScaleLogic,:) = [];
end
NMSKeypoint.Count = length(NMSKeypoint.LocationIndex);

[SphereX, SphereY, SphereZ] = sphere(100);
SphereX = reshape(SphereX,[],1);
SphereY = reshape(SphereY,[],1);
SphereZ = reshape(SphereZ,[],1);

[Value, LocationIndex] = sort(NMSKeypoint.Scale,'descend');

for i = 1 : 500 %NMSKeypoint.Count
	SphereAtPoint = bsxfun(@plus, NMSKeypoint.Scale(LocationIndex(i))*[SphereX, SphereY, SphereZ], NMSKeypoint.Location(LocationIndex(i),:));
    hh = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
    set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
global ProjectRoot;


ModelFolder = 'armadillo/';
Model = 'Armadillo_e1_100000';


FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud = findMeshNormals(PointCloud);
FileLocationNeighbors = strcat('/media/andrew/WDRhodes/diffusiondata/',ModelFolder,'neighbors/');
load(strcat(FileLocationNeighbors, Model,'_Neighbors.mat'), 'Neighbors')
PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);






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
clear
global ProjectRoot;

ModelFolder = 'buddha/';
Model = 'Buddha_e1_50000';


FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud = findMeshNormals(PointCloud);
PointCloud = findMeshResolution(PointCloud, 'Model');
FileLocationNeighbors = strcat('/media/andrew/WDRhodes/diffusiondata/',ModelFolder,'neighbors/');
load(strcat(FileLocationNeighbors, Model,'_Neighbors.mat'), 'Neighbors')
PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);


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
clear
global ProjectRoot;
ModelFolder = 'dragon/';
Model = 'Dragon_e1_50000';


FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud = findMeshNormals(PointCloud);
FileLocationNeighbors = strcat('/media/andrew/WDRhodes/diffusiondata/',ModelFolder,'neighbors/');
load(strcat(FileLocationNeighbors, Model,'_Neighbors.mat'), 'Neighbors')
PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);


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
clear
global ProjectRoot;

ModelFolder = 'itokawa/';
Model = 'Itokawa_e1_80000';


FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud = findMeshNormals(PointCloud);
FileLocationNeighbors = strcat('/media/andrew/WDRhodes/diffusiondata/',ModelFolder,'neighbors/');
load(strcat(FileLocationNeighbors, Model,'_Neighbors.mat'), 'Neighbors')
PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);



sunvector = [0;-1;-1];
sunvector = sunvector / norm(sunvector);
intensity = 0.9*PointCloud.Normal*sunvector;
intensity(intensity<0) = 0;
intensity(isnan(intensity)) = 0;


FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'SignalNoise/Mesh_te8/');
load(strcat(FileLocation,'NMSKeypoint.mat'),'NMSKeypoint')



figure
colormap('gray')
% trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), 'EdgeColor', 'none');
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), intensity, 'EdgeColor', 'none');
view(0, -40)
axis equal
axis off
hold on















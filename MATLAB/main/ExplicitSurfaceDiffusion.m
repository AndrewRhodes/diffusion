% Andrew Rhodes
% ASEL
% October 2017


close all
clear
clc




%% Additional Paths

addpath('~/GitProjects/pose/MATLAB_PointCloudDescriptors/OURCVFH/models/')
addpath('~/GitProjects/matlab-utilities/')
addpath('~/Desktop/MLIDAR-1.0/MATLAB_Modules/')
addpath('~/AFOSR/Ashish/CS3 Code/')


%% Load Model

[vertices, faces] = read_ply('CleanItokawa.ply');
Model = pcread('CleanItokawa.ply');

PointCloud.Location = vertices;
PointCloud.Face = faces;
PointCloud.Normal = Model.Normal;
PointCloud.Count = length(vertices);
PointCloud.FaceCount = size(faces, 1);
PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceNormal = findFaceNormals(PointCloud.Location, PointCloud.Face);



%% Downsample the Mesh

PC.faces=PointCloud.Face;
PC.vertices = PointCloud.Location;

PC2 = reducepatch(PC, 90000);

clear PointCloud

PointCloud.Face = PC2.faces;
PointCloud.Location = PC2.vertices;

[VerticesOut, FacesOut] = clearMeshDuplicates(PointCloud.Location, PointCloud.Face );

clear PointCloud


PointCloud.Face = FacesOut;
PointCloud.Location = VerticesOut;
PointCloud = findMeshNormals(PointCloud);
PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.Count = length(vertices);
PointCloud.FaceCount = size(PointCloud.Face, 1);


% save_off(PointCloud.Location, PointCloud.Face, 'Itokawa90000.off')


%% Load Laplace Beltrami Operator
load('AItokawa.mat')
load('LapMatMeshWeightsItokawa.mat')
load('hItokawa.mat')

% load('AKepler90000.mat')
% load('LapMatMeshWeightsKepler90000.mat')
% load('hKepler90000.mat')

dt = (h/2)^2;


Alength = length(A);

A1 = sparse(1:Alength,1:Alength, 1./A);


LBM = A1 * LapMatMeshWeights;






%% Estimate Guassian Scale Paramter

% Number of diffusion steps
MaxLevel = 30;

% maxsample = max(eigs(-LBM));
maxsample = 5.1032;

ws = 0 : 0.001 : maxsample;

NumSample = length(ws);

ScaleParameterFrequency = zeros(MaxLevel,1);
ScaleParameterFrequency(1,1) = inf;
ScaleParamterSpatial = zeros(MaxLevel,1);
CutoffFrequency = zeros(MaxLevel,1);
 

% Cutoff Frequency scalar
cut = sqrt(log(2));

% Decibal cutoff
db3 = 1/sqrt(2);

%Compute freq squared
ws2 = (ws.^2)';
ws4 = (ws.^4)';

H = zeros(NumSample, MaxLevel);
H(:,1) = 1./ (ones(NumSample,1) + 1 * ws2);
h = zeros(NumSample, 1);

% Get transfer functions
for i = 2 : MaxLevel+1
    
    % Transfer function
    h(:,1) = 1./ (ones(NumSample,1) + 1 * ws2);
    H(:,i) = h(:,1) .* H(:,i-1);
    
    % Find the frequency at the cutoff values
    CutoffFrequency(i,1) = interp1(H(:,i), ws, db3);
    
    % Change cutoff frequency to scale parameter
    ScaleParameterFrequency(i,1) = CutoffFrequency(i,1) / cut;
    
    ScaleParamterSpatial(i,1) = sqrt(dt) / ScaleParameterFrequency(i,1);
end

ScaleParamterSpatial(1) = [];


%% Setup Surface Normals

SurfaceNormals = zeros(PointCloud.LocationCount, 3, MaxLevel);
SurfaceNormals(:,:,1) = PointCloud.Normal;



%% Precompute Operator Inverse

SparseIdentity = speye(PointCloud.LocationCount, PointCloud.LocationCount);

ILBM = (SparseIdentity - dt * LBM);



%% Implicit Euler Smoothing Scheme Normals


% clearvars -except ILBM MaxLevel SurfaceNormals

WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, MaxLevel-1));

for i = 1 : MaxLevel
    
    [a(:,1), flag] = bicg(ILBM, SurfaceNormals(:,1,i), 1e-10, 50);
    [a(:,2), flag] = bicg(ILBM, SurfaceNormals(:,2,i), 1e-10, 50);
    [a(:,3), flag] = bicg(ILBM, SurfaceNormals(:,3,i), 1e-10, 50);

    SurfaceNormals(:,1:3,i+1) = bsxfun(@rdivide, a, sqrt(sum(a.^2,2)) );
    
    waitbar(i/MaxLevel, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, MaxLevel-1));
    
end

waitbar(i/MaxLevel, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)

% 
% save SurfaceNormalsKepler SurfaceNormals

%% Setup Curvature Values

Neighbors = findAdjacentNeighbors(PointCloud.Face, PointCloud.Location);

MoreNeighbors = cell(PointCloud.LocationCount,1);

for i = 1 : PointCloud.LocationCount
    CurrentNeighbors = unique([Neighbors{Neighbors{i,1}}]);
%     CurrentNeighbors = unique([Neighbors{[Neighbors{[Neighbors{Neighbors{i,1}}]}]}]);
    CurrentNeighbors(CurrentNeighbors == i) = [];
    
    MoreNeighbors{i,1} = CurrentNeighbors;
end


[PK1, PK2, PD1, PD2, MK, GK] = findPointCurvatures(PointCloud.Location, PointCloud.Normal, Neighbors);

SurfaceCurvatures = zeros(PointCloud.LocationCount, 1, MaxLevel);
% SurfaceCurvatures(:,:,1:30) = SurfaceCurvaturesA;
SurfaceCurvatures(:,:,1) = MK;


%% Implicit Euler Smoothing Scheme Curvature


WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, MaxLevel-1));

for i = 1 : MaxLevel
    
    SurfaceCurvatures(:,1,i+1) = bicg(ILBM, SurfaceCurvatures(:,:,i), 1e-10, 50);
    
    waitbar(i/MaxLevel, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, MaxLevel-1));
    
end

waitbar(i/MaxLevel, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)

% save SurfaceCurvaturesKepler SurfaceCurvatures

%% Estimated Laplacians (Difference of Gaussian DoG)
% Defined by Fadaifard
% Scale Normalized Laplacian, and scale Invariant Laplacian

LaplacianScale = zeros(PointCloud.LocationCount, 1, MaxLevel-1);
LaplacianScaleNormalized = zeros(PointCloud.LocationCount, 1, MaxLevel-1);
LaplaceScaleInvariant = zeros(PointCloud.LocationCount, 1, MaxLevel-1);


for i = 1 : MaxLevel-1
    
    % Estimated Laplacian
    LaplacianScale(:,1,i) = ( SurfaceCurvatures(:,1,i+1) - SurfaceCurvatures(:,1,i) ) / (ScaleParamterSpatial(i+1,1)^2-ScaleParamterSpatial(i,1)^2);
    
    % Scale Normalized Laplacian
    LaplacianScaleNormalized(:,1,i) = 2 * ( SurfaceCurvatures(:,1,i+1) - SurfaceCurvatures(:,1,i) ) * (ScaleParamterSpatial(i,1)^2 / (ScaleParamterSpatial(i+1,1)^2-ScaleParamterSpatial(i,1)^2) );
    
    % Scale Invariant Laplace
    Fbar = (sum(SurfaceCurvatures(:,1,i)) / PointCloud.LocationCount) * ones(PointCloud.LocationCount, 1) ;
    
    SigmaL = sqrt(sum((LaplacianScale(:,1,i)-Fbar).^2)) / sqrt(PointCloud.LocationCount);
    %  
    LaplaceScaleInvariant(:,1,i) = ( LaplacianScale(:,1,i) - Fbar ) / SigmaL;
    
end



%% Extrema Detection

% Use this laplacian descriptor
LaplaceDOG = LaplacianScaleNormalized;


NumFeatures = 0;

for i = 1 : PointCloud.LocationCount
    
    CurrentNeighbors = Neighbors{i,1};
    
    for j = 2 : MaxLevel - 2
        
        CurrentValue = LaplaceDOG(i,1,j);
        
        SurroundingValues = [LaplaceDOG([i,CurrentNeighbors],1,j-1);
            LaplaceDOG(CurrentNeighbors,1,j);
            LaplaceDOG([i,CurrentNeighbors],1,j+1)];
        
        % Check Both if it is maximum or minimum
        if all(CurrentValue > SurroundingValues)
            NumFeatures = NumFeatures + 1;
            % Maximum
            FeaturePoint.Scale(NumFeatures, 1) = ScaleParamterSpatial(j);
            FeaturePoint.Location(NumFeatures, 1) = i;
        end
        if all(CurrentValue < SurroundingValues)
            % Minimum
             NumFeatures = NumFeatures + 1;
            FeaturePoint.Scale(NumFeatures, 1) = ScaleParamterSpatial(j);
            FeaturePoint.Location(NumFeatures, 1) = i;
        end
        
    end
    
end


%% Find Features of Max Scale

[SortValue, SortOrder] = sort(FeaturePoint.Scale,'descend');

UniqueSortValue = flipud(unique(SortValue));

PointSizes = 2*round(10*unique(SortValue));

FeaturePointSizes = zeros(length(SortValue),1);

colors = ['r','g','b','y','c','m','w','r','g','b','y','c','m','w'];

for i = 1 : 50
    
    FeaturePointSizes(i,1) = PointSizes(SortValue(i) == UniqueSortValue);
    FeaturePointColor(i,1) = colors(SortValue(i) == UniqueSortValue);
    
end




%% Plot Features

[SphereX, SphereY, SphereZ] = sphere(100);

SphereX = reshape(SphereX,[],1);
SphereY = reshape(SphereY,[],1);
SphereZ = reshape(SphereZ,[],1);


sunvector = [0.1;-1;-1];
sunvector = sunvector / norm(sunvector);

% FaceNormalsInCamera = SurfaceNormals(:,:,1);

FaceNormalsInCamera = PointCloud2.Normal;


intensity = 0.9*FaceNormalsInCamera*sunvector;
intensity(intensity<0) = 0;


PlotFeatures = 1:50;


figure
colormap('gray')
% trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), intensity, 'EdgeColor', 'none');
trisurf(PointCloud2.Face, PointCloud2.Location(:,1), PointCloud2.Location(:,2), PointCloud2.Location(:,3), intensity, 'EdgeColor', 'none');
% view(-180, 0)
view(15, 5)
axis equal
axis off
hold on



for i = 1 : length(PlotFeatures)
    
    SphereAtPoint = bsxfun(@plus, SortValue(PlotFeatures(i))*[SphereX, SphereY, SphereZ], PointCloud.Location(FeaturePoint.Location(SortOrder(PlotFeatures(i))), :));

    hh = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
%     set(hh, 'FaceColor',FeaturePointColor(PlotFeatures(i)), 'EdgeColor','none', 'FaceAlpha',0.5)
    set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)

    
    
%     pause
end






%% Create Scale Space Representation



Vector = -PointCloud.Normal(FeaturePoint.Location(SortOrder(1)),:);
Axis = cross(Vector, [0,0,1]);
e_theta = Axis/norm(Axis);
theta = asin(norm(Axis));

Transform = v_to_T(theta*e_theta);

% (Transform*PointCloud.Normal([MoreNeighbors{FeaturePoint.Location(SortOrder(1))},FeaturePoint.Location(SortOrder(1))],:)')'


a=bsxfun(@minus, (Transform*PointCloud.Location([MoreNeighbors{FeaturePoint.Location(SortOrder(1))},FeaturePoint.Location(SortOrder(1))],:)')', (Transform*PointCloud.Location(FeaturePoint.Location(SortOrder(1)),:)')')

imshow(a(:,1:2))

for level = 1 : MaxLevel

    
    PointCloud.Location(FeaturePoint.Location(SortOrder(1)),:)
    PointCloud.Normal(FeaturePoint.Location(SortOrder(1)),:)
    
    CurrentLevelCurvatures = SurfaceCurvatures(MoreNeighbors{FeaturePoint.Location(SortOrder(1))}, 1, level);
    
    hist(CurrentLevelCurvatures)
    xlim([2e-3, 4e-3])
    pause
    
    
    PointCloud.Location(MoreNeighbors{FeaturePoint.Location(SortOrder(PlotFeatures(i)))}, 1),...
            PointCloud.Location(MoreNeighbors{FeaturePoint.Location(SortOrder(PlotFeatures(i)))}, 2),...
            PointCloud.Location(MoreNeighbors{FeaturePoint.Location(SortOrder(PlotFeatures(i)))}, 3)

end


























% Andrew Rhodes
% October 2017
% Diffusion on a mesh using Colin B. Macdonald's closest point algorithm,
% Bertalmio's level set approach, and my implicit Euler method



close all
clear
clc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the triangular surface 
Model = 'CleanItokawa.ply';
% Model = 'MGS_launch.3ds.ply';
% Model = 'MGS_moi.ply';

% Spacing of implicit surface points
spacing = 0.5; 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Defined variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% degree of interpolation
porder = 3; 
% dimension of the surface
dim = 3;

% Order of numerical differncing. 1 for 2nd order, 2 for 4th order
Lorder = 2;

% define bandwidth around surface
bandwidth = 1.00001*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));


tauImplicit = spacing / 4;
MaxTauImplicit = 20;
NumStepsImplicit = MaxTauImplicit / tauImplicit;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional Paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('~/AFOSR/Ashish/CS3 Code')
addpath('~/GitProjects/pose/MATLAB_PointCloudDescriptors/OURCVFH/models/')
addpath('~/GitProjects/matlab-utilities/')
addpath('~/GitProjects/mesh_resampling_toolbox/MATLAB_Modules/')
addpath(genpath('~/Documents/Software/cp_matrices/'))
addpath('~/Desktop/MLIDAR-1.0/MATLAB_Modules/')
addpath('~/Downloads/')
addpath(genpath('~/GitProjects/pose/MLIDAR-2.2/'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[vertices, faces] = read_ply(Model);
PCModel = pcread(Model);

PointCloud.Location = vertices;
PointCloud.Face = faces;
% PointCloud.Normal = PCModel.Normal;
% PointCloud.Count = length(vertices);
% PointCloud.FaceCount = size(faces, 1);
% PointCloud.LocationCount = size(PointCloud.Location,1);


% PC.faces=PointCloud.Face;
% PC.vertices = PointCloud.Location;
% 
% PC2 = reducepatch(PC, 40000);
% 
% clear PointCloud
% 
% PointCloud.Face = PC2.faces;
% PointCloud.Location = PC2.vertices;
% 
% [vertices, faces, vertexNormals, faceNormals, adjFaceIndices, octreeResolution] = mesh_resample(Model);

[VerticesOut, FacesOut] = clearMeshDuplicates(PointCloud.Location, PointCloud.Face );
% 
% clear PointCloud
% 
% 
PointCloud.Face = FacesOut;
PointCloud.Location = VerticesOut;
PointCloud = findMeshNormals(PointCloud);
PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.Count = length(PointCloud.Location);
PointCloud.FaceCount = size(PointCloud.Face, 1);

PointCloud.FaceNormal = findFaceNormals(PointCloud.Location, PointCloud.Face);
PointCloud.FaceArea = findFaceAreas(PointCloud.Location, PointCloud.Face);


clear vertices faces Model PCModel PC PC2 VerticesOut FacesOut
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make starting verticies for implicit grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MinPoint = round(min(PointCloud.Location) - bandwidth - spacing, 1);
MaxPoint = round(max(PointCloud.Location) + bandwidth + spacing, 1);
% Find Closet Point

% XYZ = MinPoint - spacing + IJK * spacing

% CP is closest point on triangular surface, not necessarily a vertex,
% could exist on a face.

[IJK,DIST,CP,XYZ,CPFACE] = tri2cp(PointCloud.Face, PointCloud.Location, spacing, MinPoint, porder, 1);


NumGridPoints = length(IJK);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expand Surface Data To Implicit Space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Neighbors = findAdjacentNeighbors(PointCloud.Face, PointCloud.Location);

[PK1, PK2, PD1, PD2, MK, GK] = findPointCurvatures(PointCloud.Location, PointCloud.Normal, Neighbors);

PointCloud.Signal = MK;

x1d = (MinPoint(1):spacing:MaxPoint(1))';
y1d = (MinPoint(2):spacing:MaxPoint(2))';
z1d = (MinPoint(3):spacing:MaxPoint(3))';

BandSearchSize = [length(y1d), length(x1d), length(z1d)];

Band = sub2ind(BandSearchSize, IJK(:,2), IJK(:,1), IJK(:,3));


% Interpolate vertex values to closest points on surface
FaceInterpolateWeights = interpBarycenterTriangle(PointCloud, CP, CPFACE);


Signal = FaceInterpolateWeights * PointCloud.Signal;




% Create L, E, M
L = laplacian_3d_matrix(x1d,y1d,z1d, Lorder, Band);


Eplot = interp3_matrix(x1d, y1d, z1d, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), porder, Band);
Ecp = interp3_matrix(x1d, y1d, z1d, CP(:,1), CP(:,2), CP(:,3), porder, Band);

% Eplot = interpLagrange3D(BandSearchSize, MinPoint, PointCloud.Location, porder, Band, spacing);

% Ecp = interpLagrange3D(BandSearchSize, MinPoint, CP, porder, Band, spacing);

M = lapsharp(L, Ecp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precompute the scale parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the scale parameters
maxsample = 10; 
ws = 0 : 0.01 : maxsample;
NumSample = length(ws);
ScaleParameterSpatialImplicit = zeros(NumStepsImplicit,2);
cut = sqrt(log(2));
db3 = 1/sqrt(2);
ws2 = (ws.^2)';
H = ones(NumSample,1) ;
h = 1 ./ ( ones(NumSample,1) + tauImplicit/spacing^2 * ws2 );

for i = 1 : NumStepsImplicit - 1
    
    % Transfer function
    % 1st order
    H = H .* h;
    % Find the frequency at the cutoff values
    [uH, indH] = unique(H);
    CutoffFrequencyImplicit = interp1(uH, ws(indH), db3);
    % Change cutoff frequency to scale parameter
    ScaleParameterFrequencyImplicit = CutoffFrequencyImplicit / cut;   
    ScaleParameterSpatialImplicit(i+1,1) =  spacing / ScaleParameterFrequencyImplicit;
    ScaleParameterSpatialImplicit(i+1,2) = sqrt(2*i*tauImplicit);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Diffusion on Implicit Mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Perform diffusion

SignalImplicit = zeros(PointCloud.LocationCount, NumStepsImplicit);
SignalImplicit(:,1) = PointCloud.Signal;


ItM = speye(size(M)) - tauImplicit * M;

WaitBar = waitbar(0, sprintf('Implicit Diffusion %i of %i', 0, NumStepsImplicit-1));

for i = 1 : NumStepsImplicit - 1
    

    [SignalNew, flag, relres] = bicg(ItM, Signal, 1e-10, 30);

%     [SignalNew, flag, relres] = gmres(ItM, Signal, [], 1e-10, 30);

    if flag
        flag
        relres
    end
    
    
    Signal = SignalNew;
    
    SignalImplicit(:,i+1) = Eplot * SignalNew;
    
    waitbar(i/NumStepsImplicit, WaitBar, sprintf('Implicit Diffusion %i of %i', i, NumStepsImplicit-1));

end

waitbar(i/NumStepsImplicit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimated Laplacians (Difference of Gaussian DoG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defined by Fadaifard
% Scale Normalized Laplacian, and scale Invariant Laplacian

LaplacianScale = zeros(PointCloud.LocationCount, NumStepsImplicit-1);
LaplacianScaleNormalized = zeros(PointCloud.LocationCount, NumStepsImplicit-1);
LaplaceScaleInvariant = zeros(PointCloud.LocationCount, NumStepsImplicit-1);


for i = 1 : NumStepsImplicit-1
    
    % Estimated Laplacian
    LaplacianScale(:,i) = ( SignalImplicit(:,i+1) - SignalImplicit(:,i) ) / (ScaleParameterSpatialImplicit(i+1,1)^2-ScaleParameterSpatialImplicit(i,1)^2);
    
    % Scale Normalized Laplacian
    LaplacianScaleNormalized(:,i) = 2 * ( SignalImplicit(:,i+1) - SignalImplicit(:,i) ) * (ScaleParameterSpatialImplicit(i,1)^2 / (ScaleParameterSpatialImplicit(i+1,1)^2 - ScaleParameterSpatialImplicit(i,1)^2) );
    
    % Scale Invariant Laplace
    Fbar = (sum(SignalImplicit(:,i)) / PointCloud.LocationCount) * ones(PointCloud.LocationCount, 1) ;
    
    SigmaL = sqrt(sum((LaplacianScale(:,i)-Fbar).^2)) / sqrt(PointCloud.LocationCount);
    %  
    LaplaceScaleInvariant(:,i) = ( LaplacianScale(:,i) - Fbar ) / SigmaL;
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extrema Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this laplacian descriptor
LaplaceDOG = LaplacianScaleNormalized;

FeaturePoint.Scale = zeros(11000,2);
FeaturePoint.Location = zeros(11000,1);

NumFeatures = 0;

for i = 1 : PointCloud.LocationCount
    
    CurrentNeighbors = Neighbors{i,1};
    
%     for j = 2 : 150 - 2
    for j = 2 : NumStepsImplicit - 2
        
        CurrentValue = LaplaceDOG(i,j);
        
        SurroundingValues = [LaplaceDOG([i,CurrentNeighbors],j-1);
            LaplaceDOG(CurrentNeighbors,j);
            LaplaceDOG([i,CurrentNeighbors],j+1)];
        
        % Check Both if it is maximum or minimum
        if all(CurrentValue > SurroundingValues)
            NumFeatures = NumFeatures + 1;
            % Maximum
            FeaturePoint.Scale(NumFeatures, :) = ScaleParameterSpatialImplicit(j,:);
            FeaturePoint.Location(NumFeatures, 1) = i;
        end
        if all(CurrentValue < SurroundingValues)
            % Minimum
             NumFeatures = NumFeatures + 1;
            FeaturePoint.Scale(NumFeatures, :) = ScaleParameterSpatialImplicit(j,:);
            FeaturePoint.Location(NumFeatures, 1) = i;
        end
        
    end
    
end

FeaturePoint.Scale(FeaturePoint.Location == 0,:) = [];
FeaturePoint.Location(FeaturePoint.Location == 0) = [];

% save FeaturePointMGSMOI FeaturePoint




%% Find Features of Max Scale

[SortValue, SortOrder] = sort(FeaturePoint.Scale(:,1),'descend');


MediumFeatures = find(SortValue>1 & SortValue<4);

%% Plot Features

[SphereX, SphereY, SphereZ] = sphere(100);

SphereX = reshape(SphereX,[],1);
SphereY = reshape(SphereY,[],1);
SphereZ = reshape(SphereZ,[],1);


% sunvector = [0.1;.1;1];
sunvector = [1;1;0];
sunvector = sunvector / norm(sunvector);

% FaceNormalsInCamera = SurfaceNormals(:,:,1);

FaceNormalsInCamera = PointCloud.Normal;


intensity = 0.9*FaceNormalsInCamera*sunvector;
intensity(intensity<0) = 0;


PlotFeatures = 11000:11100;

% PlotFeatures = 10070:10080;

PlotOrder = [3,1,2];

figure
colormap('gray')
trisurf(PointCloud.Face, PointCloud.Location(:,PlotOrder(1)), PointCloud.Location(:,PlotOrder(2)), PointCloud.Location(:,PlotOrder(3)), intensity, 'EdgeColor', 'none');
% view(-180, 0)
view(15, 5)
axis equal
axis off
hold on
% xlabel('x')
% ylabel('y')
% zlabel('z')
view(-187,-61)


for i = PlotFeatures%1 : length(PlotFeatures)
    
%     SphereAtPoint = bsxfun(@plus, SortValue(PlotFeatures(i))*[SphereX, SphereY, SphereZ], PointCloud.Location(FeaturePoint.Location(SortOrder(PlotFeatures(i))), :));
      SphereAtPoint = bsxfun(@plus, SortValue(i)*[SphereX, SphereY, SphereZ], PointCloud.Location(FeaturePoint.Location(SortOrder(i)), :));

%     SphereAtPoint = bsxfun(@plus, SortValue(MediumFeatures(i))*[SphereX, SphereY, SphereZ], PointCloud.Location(FeaturePoint.Location(SortOrder(MediumFeatures(i))), :));
%     SphereAtPoint = bsxfun(@plus, [SphereX, SphereY, SphereZ], PointCloud.Location(FeaturePoint.Location(SortOrder(MediumFeatures(i))), :));

    hh = surfl(reshape(SphereAtPoint(:,PlotOrder(1)),101,101), reshape(SphereAtPoint(:,PlotOrder(2)),101,101), reshape(SphereAtPoint(:,PlotOrder(3)),101,101));
%     set(hh, 'FaceColor',FeaturePointColor(PlotFeatures(i)), 'EdgeColor','none', 'FaceAlpha',0.5)
    set(hh, 'FaceColor','b', 'EdgeColor','none', 'FaceAlpha',0.5)
%     view(170,70)
%     view(136,13)
%     view(34, 19)
%     drawnow
    
%     pause
end



for i = 1 : length(PlotFeatures)
    
    SphereAtPoint = bsxfun(@plus, SortValue(end + 1 - PlotFeatures(i))*[SphereX, SphereY, SphereZ], PointCloud.Location(FeaturePoint.Location(SortOrder(end + 1 - PlotFeatures(i))), :));

    hh = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
%     set(hh, 'FaceColor',FeaturePointColor(PlotFeatures(i)), 'EdgeColor','none', 'FaceAlpha',0.5)
    set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',1)

    
    
%     pause
end



















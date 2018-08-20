% Andrew Rhodes
% ASEL
% March 2018


close all
clear

clc

global ProjectRoot; % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


alpha = 1;

options.rho = 4;
options.dtype = 'geodesic';
% options.dtype = 'euclidean';

ShowPlot = 0;
Model = 'bunny/Bunny_e1';

BDF = 2;
tauFraction = 1/10;
tauNumerator = 250;
DoGNormalize = 'DoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'



NumIter = 50;

t_scale = 0.7;
t_DoG = 0.9;
t_range = 3;

NoiseVec = [0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model File Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(Model,'.ply');
FileNameModelOff = strcat(Model,'.off');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

% PointCloud.Location = PointCloud.Location.*681.65;

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
PointCloud = findMeshResolution(PointCloud, 'Model');
PointCloud = findMeshNormals(PointCloud)
NormalRotations = findNormalsRotation(PointCloud.Normal);

% ply_write('Bunny_e1', PointCloud.Face, PointCloud.Location)
% save_off(PointCloud.Location, PointCloud.Face, 'Bunny_e1.off')


% % % % % % % % % %
tau = PointCloud.Resolution * tauFraction;
MaxTau = tauNumerator / PointCloud.Resolution;
NumSteps = round(MaxTau);
% % % % % % % % % %

load('BunnyCurvature_e1.mat')
MK = Curvature;

load('Bunny_e1_Neighbors.mat')
PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);

%[Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);

% [PK1, PK2, PD1, PD2, MK, GK] = findPointCurvatures(PointCloud, NormalRotations, Neighbors.Connect);
clear PK1 PK2 PD1 PD2 PK2 GK NeighborFaces NormalRotations
stdMK = std(MK);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ItL = makeExplicitLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), options, BDF, tau, alpha);

load('BunnyItL_e1.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ScaleParameter = findScaleParamter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusion of Mean Curvature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : length(NoiseVec)

    for i = 27:50
        i
        PointCloud.Signal = MK + NoiseVec(j)*stdMK*rand(PointCloud.LocationCount,1);
        
        
        Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find Difference of Gaussian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Detect Extrema
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Keypoint = findKeypoint(DoG, PointCloud, ScaleParameter, Neighbors.Distance, KeypointMethod, CompareMethod);
        NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);
        
        %     SubKeypoint = findSubKeypoint(Keypoint, ScaleParameterAbsolute, DoG, PointCloud, Neighbors.Connect, NeighborFaces.Connect);
        
        
        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/bunny/Std_',num2str(NoiseVec(j)));
       FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
%       FileName = strcat('Keypoint','.mat');
        
        save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
        
       FileName = strcat('NMSKeypoint','_Iter',num2str(i),'.mat');
%        FileName = strcat('NMSKeypoint','.mat');
       
        save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
        
    end
end




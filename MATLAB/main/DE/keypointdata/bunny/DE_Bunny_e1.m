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

ShowPlot = 1;
Model = 'bunny/Bunny_e1';

BDF = 4;
tauFraction = 1/10;
DoGNormalize = 1;
tauNumerator = 250;

NumIter = 50;

t_scale = 0.7;
t_DoG = 0.9;

NoiseVec = [0.1, 0.2, 0.3];

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
PointCloud = findMeshNormals(PointCloud);
NormalRotations = findNormalsRotation(PointCloud.Normal);

% ply_write('Bunny_e1', PointCloud.Face, PointCloud.Location)
% save_off(PointCloud.Location, PointCloud.Face, 'Bunny_e1.off')


% % % % % % % % % %
tau = PointCloud.Resolution * tauFraction;
MaxTau = tauNumerator / PointCloud.Resolution;
NumSteps = round(MaxTau);
% % % % % % % % % %


[Neighbors, NeighborFaces] = findAdjacentNeighbors(PointCloud);

[PK1, PK2, PD1, PD2, MK, GK] = findPointCurvatures(PointCloud, NormalRotations, Neighbors.Connect);
clear PK1 PK2 PD1 PD2 PK2 GK NeighborFaces NormalRotations
stdMK = std(MK);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ItL = makeExplicitLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), options, BDF, tau, alpha);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ScaleParameter = findScaleParamter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusion of Mean Curvature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : length(NoiseVec)
    if NoiseVec(j) == 0.3
        Iterand = [1:NumIter/2];
    else
        Iterand = [1:NumIter];
    end
    for i = Iterand
        
        PointCloud.Signal = MK + NoiseVec(j)*stdMK*rand(PointCloud.LocationCount,1);
        
        
        Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find Difference of Gaussian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Detect Extrema
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Keypoint = findKeypoint(DoG, ScaleParameter, Neighbors.Distance, 'Old');

        NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_DoG);
        
        %     SubKeypoint = findSubKeypoint(Keypoint, ScaleParameterAbsolute, DoG, PointCloud, Neighbors.Connect, NeighborFaces.Connect);
        
        
        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/bunny/Std_',num2str(NoiseVec(j)),'_new/');
        FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
        %     FileName = strcat('Keypoint','.mat');
        
        save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
        
        FileName = strcat('NMSKeypoint','_Iter',num2str(i),'.mat');
        %     FileName = strcat('NMSKeypoint','.mat');
        
        save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
        
    end
end






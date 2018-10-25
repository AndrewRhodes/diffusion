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


alpha = 1/2;

options.rho = 6;
options.dtype = 'geodesic';
% options.dtype = 'euclidean';
options.htype = 'psp';


ModelFolder = 'itokawa/';
Model = 'Itokawa_e1_80000';


BDF = 2;
tauFraction = 1/10;
NumIter = 10;
NumSteps = 5000;

KeypointSearchMethod = 'Global'; % 'Local', 'Global'
Radius = 1;


NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model File Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');
FileNameModelOff = strcat(ModelFolder,Model,'.off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
PointCloud = findMeshResolution(PointCloud, 'Model');
PointCloud = findMeshNormals(PointCloud);

% [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);

load(strcat(Model,'_Neighbors.mat'),'Neighbors')
PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);
options.hs = PointCloud.Resolution/2;


% % % % % % % % % %
tau = options.hs^2/4;
% tau = PointCloud.Resolution * tauFraction;
% MaxTau = tauNumerator / PointCloud.Resolution;
% NumSteps = round(MaxTau);
% NumSteps = tauNumerator;
% % % % % % % % % %

load(strcat(Model,'_Curvature.mat'),'Curvature')

MK = Curvature;


% [PK1, PK2, PD1, PD2, MK, GK] = findPointCurvatures(PointCloud, NormalRotations, Neighbors.Connect);
% clear PK1 PK2 PD1 PD2 GK NeighborFaces
Quants = quantile(MK, [0.25,0.5,0.75]);
MKQuant = MK;
OutOfBounds = (MKQuant > (Quants(3) + 1.5*(Quants(3)-Quants(1)))) | (MKQuant < (Quants(1) - 1.5*(Quants(3)-Quants(1))));
MKQuant(OutOfBounds) = [];
stdMK = std(MKQuant);
% stdMK = std(MK);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(strcat(Model,'_ItL_rho',num2str(options.rho),'_dtype_',options.dtype,'.mat'),'ItL')
% ItL = makeExplicitLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), options, BDF, tau, alpha);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ScaleParameter = findScaleParamter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusion of Mean Curvature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : length(NoiseVec)
    
    for i = 1 : NumIter
        i
        PointCloud.Signal = MK + NoiseVec(j)*stdMK*randn(PointCloud.LocationCount,1);
        
        Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find Difference of Gaussian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [DoG, AbsDoG] = buildDoGNew(Signal, ScaleParameter);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Detect Extrema
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [Keypoint, Zeropoint] = findKeypointNew(PointCloud, ScaleParameter, Neighbors.Distance, DoG, AbsDoG);
        
        NewKeypoint = checkKeypointSphere(PointCloud, ScaleParameter, Keypoint, Zeropoint, Radius, KeypointSearchMethod);
        
        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'SignalNoise/NewRun/Std_',num2str(NoiseVec(j)),'/');
        FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
        
        FileName = strcat('NewKeypoint','_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'NewKeypoint', '-v7.3')
        
        FileName = strcat('Zeropoint','_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'Zeropoint', '-v7.3')

	clear Keypoint Zeropoint NewKeypoint Signal DoG AbsDoG
        
    end
end




for i = 1
    i
    PointCloud.Signal = MK;
    
    Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find Difference of Gaussian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [DoG, AbsDoG] = buildDoGNew(Signal, ScaleParameter);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detect Extrema
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Keypoint, Zeropoint] = findKeypointNew(PointCloud, ScaleParameter, Neighbors.Distance, DoG, AbsDoG);
    
    NewKeypoint = checkKeypointSphere(PointCloud, ScaleParameter, Keypoint, Zeropoint, Radius, KeypointSearchMethod);
    
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'SignalNoise/NewRun/');
    FileName = strcat('Keypoint','.mat');
    save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
    
    FileName = strcat('NewKeypoint','.mat');
    save(fullfile(FileLocation, FileName), 'NewKeypoint', '-v7.3')
    
    FileName = strcat('Zeropoint','.mat');
    save(fullfile(FileLocation, FileName), 'Zeropoint', '-v7.3')
    
    
end



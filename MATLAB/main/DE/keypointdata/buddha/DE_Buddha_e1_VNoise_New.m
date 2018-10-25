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

ModelFolder = 'buddha/';
Model = 'Buddha_e1_50000';

BDF = 2;
tauFraction = 1/10;
tauNumerator = 3000;

KeypointSearchMethod = 'Global'; % 'Local', 'Global'
Radius = 1;

NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];


for j = 1 : length(NoiseVec)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Model File Location
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    FileLocationModel = strcat(ProjectRoot,'/models/object/');
    FileNameModelPly = strcat(ModelFolder,Model,'_sigma',num2str(j),'.ply');
    FileNameModelOff = strcat(ModelFolder,Model,'_sigma',num2str(j),'.off');
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load the Model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
    
    PointCloud.LocationCount = size(PointCloud.Location,1);
    PointCloud.FaceCount = size(PointCloud.Face, 1);
    PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
    PointCloud = findMeshResolution(PointCloud, 'Model');
    PointCloud = findMeshNormals(PointCloud);
       
    
    % % % % % % % % % %
    tau = PointCloud.Resolution * tauFraction;
    MaxTau = tauNumerator / PointCloud.Resolution;
    NumSteps = round(MaxTau);
    NumSteps = tauNumerator;
    % % % % % % % % % %
    
    load(strcat(Model,'_Curvature_sigma',num2str(j),'.mat'),'Curvature')
    
    MK = Curvature;
    
    % load('Neighbors.mat')
    % [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);
    % save Armadillo_e1_100000_Neighbors Neighbors
    
    load(strcat(Model,'_Neighbors.mat'),'Neighbors')
    PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup Laplace-Beltrami
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ItL = makeExplicitLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), options, BDF, tau, alpha);
    
    save(strcat(Model,'_ItL_sigma',num2str(j),'.mat'),'ItL','-v7.3')
    
    %     load(strcat(Model,'_ItL_sigma',num2str(j),'.mat'),'ItL')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scale Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ScaleParameter = findScaleParamter(tau, alpha, NumSteps, 'Laplacian', 'Natural');
    
    ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diffusion of Mean Curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    PointCloud.Signal = MK;
    
    
    Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find Difference of Gaussian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [DoG, AbsDoG] = buildDoGNew(Signal, ScaleParameter);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detect Extrema
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Keypoint, Zeropoint] = findKeypointNew(PointCloud, ScaleParameter, Neighbors.Connect, DoG, AbsDoG);
        
	NewKeypoint = checkKeypointSphere(PointCloud, ScaleParameter, Keypoint, Zeropoint, Radius, KeypointSearchMethod);

    
    
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'VertexNoise/NewRun/Std_',num2str(NoiseVec(j)),'/');
    FileName = strcat('Keypoint','_Iter1','.mat');
    save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
    
    FileName = strcat('NewKeypoint','_Iter1','.mat');
    save(fullfile(FileLocation, FileName), 'NewKeypoint', '-v7.3')
    
    FileName = strcat('Zeropoint','_Iter1','.mat');
    save(fullfile(FileLocation, FileName), 'Zeropoint', '-v7.3')
    
end







for i = 1
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Model File Location
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    FileLocationModel = strcat(ProjectRoot,'/models/object/');
    FileNameModelPly = strcat(ModelFolder,Model,'_sigma',num2str(j),'.ply');
    FileNameModelOff = strcat(ModelFolder,Model,'_sigma',num2str(j),'.off');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load the Model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
    
    PointCloud.LocationCount = size(PointCloud.Location,1);
    PointCloud.FaceCount = size(PointCloud.Face, 1);
    PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
    PointCloud = findMeshResolution(PointCloud, 'Model');
    PointCloud = findMeshNormals(PointCloud);
    

    % % % % % % % % % %
    tau = PointCloud.Resolution * tauFraction;
    MaxTau = tauNumerator / PointCloud.Resolution;
    NumSteps = round(MaxTau);
    NumSteps = tauNumerator;
    % % % % % % % % % %
    
    load(strcat(Model,'_Curvature.mat'),'Curvature')
    
    MK = Curvature;
    
    
    load(strcat(Model,'_Neighbors.mat'),'Neighbors')
    PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);
    
    load(strcat(Model,'_ItL.mat'),'ItL') 

    ScaleParameter = findScaleParamter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

    ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);
    
    
    PointCloud.Signal = MK;
    
    Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find Difference of Gaussian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [DoG, AbsDoG] = buildDoGNew(Signal, ScaleParameter);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detect Extrema
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Keypoint, Zeropoint] = findKeypointNew(PointCloud, ScaleParameter, Neighbors.Connect, DoG, AbsDoG);
        
	NewKeypoint = checkKeypointSphere(PointCloud, ScaleParameter, Keypoint, Zeropoint, Radius, KeypointSearchMethod);

    
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'VertexNoise/NewRun/');
    FileName = strcat('Keypoint','.mat');
    save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
    
    FileName = strcat('NewKeypoint','.mat');
    save(fullfile(FileLocation, FileName), 'NewKeypoint', '-v7.3')
    
    FileName = strcat('Zeropoint''.mat');
    save(fullfile(FileLocation, FileName), 'Zeropoint', '-v7.3')
    
   
end





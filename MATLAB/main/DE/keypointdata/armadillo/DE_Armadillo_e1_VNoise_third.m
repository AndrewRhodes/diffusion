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
Model = 'armadillo/Armadillo_e1_100000';

BDF = 2;
tauFraction = 1/10;
NumIter = 1;
tauNumerator = 3000;
DoGNormalize = 'DoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

t_scale = 0.7;
t_DoG = 0.9;
t_range = 3;

NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];


for i = 1
    
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
    
    PointCloud.LocationCount = size(PointCloud.Location,1);
    PointCloud.FaceCount = size(PointCloud.Face, 1);
    PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
    PointCloud = findMeshResolution(PointCloud, 'Model');
    PointCloud = findMeshNormals(PointCloud)
    

    % % % % % % % % % %
    tau = PointCloud.Resolution * tauFraction;
    MaxTau = tauNumerator / PointCloud.Resolution;
    NumSteps = round(MaxTau);
    NumSteps = tauNumerator;
    % % % % % % % % % %
    
    load('ArmadilloCurvature_e1_100000.mat')
    
    MK = Curvature;
    
    
    load('Armadillo_e1_100000_Neighbors.mat')
    PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);
    
    load('ArmadilloItL_e1_100000.mat') 

    ScaleParameter = findScaleParamter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

    ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);
    
    
    PointCloud.Signal = MK;
    
    Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find Difference of Gaussian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detect Extrema
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Keypoint = findKeypoint(DoG, PointCloud, ScaleParameter, Neighbors.Connect, KeypointMethod, CompareMethod);
    
    NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);

    
    for j = 1 : length(NoiseVec)
        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/armadillo/VertexNoise/LongRun/Std_',num2str(NoiseVec(j)));
        FileName = strcat('Keypoint','.mat');
        
        save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
        
        FileName = strcat('NMSKeypoint','.mat');
        save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
    end
end



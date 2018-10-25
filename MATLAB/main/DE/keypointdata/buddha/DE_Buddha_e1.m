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
options.htype = 'psp';
% options.dtype = 'euclidean';

ModelFolder = 'buddha/';
Model = 'Buddha_e1_50000';

BDF = 2;
% tauFraction = 1/10;
NumIter = 20;
NumSteps = 2000;
DoGNormalize = 'NLoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'




NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];

t_scale = 0.7;
t_DoG = 0.9;
t_range = 1/2;

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

% % PointCloud.Location = PointCloud.Location.*(1/PointCloud.Resolution);

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
PointCloud = findMeshResolution(PointCloud, 'Model');
PointCloud = findMeshNormals(PointCloud);

% [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);
% save Buddha_e1_50000_Neighbors Neighbors

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
% clear PK1 PK2 PD1 PD2 GK NeighborFaces NormalRotations
Quants = quantile(MK, [0.25,0.5,0.75]);
MKQuant = MK;
OutOfBounds = (MKQuant > (Quants(3) + 1.5*(Quants(3)-Quants(1)))) | (MKQuant < (Quants(1) - 1.5*(Quants(3)-Quants(1))));
MKQuant(OutOfBounds) = [];
stdMK = std(MKQuant);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileNameItL = strcat(Model,'_ItL_rho',num2str(options.rho),'_dtype_',options.dtype,'.mat');
if ~exist(FileNameItL, 'file')
    ItL = makeCotangentLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), BDF, tau, alpha);
    save(FileNameItL, 'ItL', '-v7.3')
else
    load( FileNameItL, 'ItL')
end

ItL(2,:) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ScaleParameter = findScaleParameter(sqrt(tau), alpha, NumSteps, 'Gaussian', 'Natural');
ScaleParameter = findScaleParameter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusion of Mean Curvature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : length(NoiseVec)
    
    for i = 1 : NumIter
        sprintf('Std %0.1f : %d',NoiseVec(j),i)
        PointCloud.Signal = MK + NoiseVec(j)*stdMK*randn(PointCloud.LocationCount,1);
        
        
        Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
%         Signal = performBDFDiffusion_cpp(ItL, PointCloud.Signal, NumSteps);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find Difference of Gaussian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Detect Extrema
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         Keypoint = findKeypoint(DoG, PointCloud, ScaleParameter, Neighbors.Distance, KeypointMethod, CompareMethod);
        [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Distance);
        Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
        Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
        Keypoint.Scale = ScaleParameter(Keypoint.Level);
        
        NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);
        
        
        % SubKeypoint = findSubKeypoint(Keypoint, ScaleParameterAbsolute, DoG, PointCloud, Neighbors, NeighborFaces);
        
        
        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'SignalNoise/Mesh1/Std_',num2str(NoiseVec(j)),'/');
        FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
        %              FileName = strcat('Keypoint','.mat');
        
        save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
        
        
        FileName = strcat('NMSKeypoint','_Iter',num2str(i),'.mat');
        %              FileName = strcat('NMSKeypoint','.mat');
        
        save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
        
        
    end
    
end



for i = 1
    i
    
    PointCloud.Signal = MK;
    
    Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
%     Signal = performBDFDiffusion_cpp(ItL, PointCloud.Signal, NumSteps);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find Difference of Gaussian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detect Extrema
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     Keypoint = findKeypoint(DoG, PointCloud, ScaleParameter, Neighbors.Distance, KeypointMethod, CompareMethod);
    [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Distance);
    Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
    Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
    Keypoint.Scale = ScaleParameter(Keypoint.Level);
    
    NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);
    
    %     SubKeypoint = findSubKeypoint(Keypoint, ScaleParameterAbsolute, DoG, PointCloud, Neighbors.Connect, NeighborFaces.Connect);
    
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'SignalNoise/Mesh1/');
    FileName = strcat('Keypoint','.mat');
    
    save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
    
    FileName = strcat('NMSKeypoint','.mat');
    save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
    
end














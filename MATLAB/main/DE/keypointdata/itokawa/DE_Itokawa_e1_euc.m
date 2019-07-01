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
options.dtype = 'euclidean'; % 'euclidean', 'geodesic' %
options.htype = 'ddr'; % 'psp', 'ddr'

Destination = 'Mesh_rho4_ddr_euc' %'Mesh_rho4_ddr_geo_t0.25'
ModelFolder = 'itokawa/';
Model = 'Itokawa_e1_80000';

NumIter = 15;
DoGNormalize = 'DoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

l_range1 = 1/2;
l_range2 = 2;
k = 1.2;%2^(1/4);
l_scale = 1/k;
l_ebar = 80;


set_ltau = @(e_bar) 4*e_bar;
set_hs = @(e_bar) 2*e_bar^(1/5);
set_t0 = @(e_bar) e_bar/4 ;
set_NumSteps = @(e_bar, t_0) ceil(log((l_ebar*e_bar)^2 / (2*alpha*t_0)) / (2*log(k)));


NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model File Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');
FileNameModelOff = strcat(ModelFolder,Model,'.off');
FileLocationWD = '/media/andrew/WDRhodes/diffusiondata/';

TmpLocation = strcat(ProjectRoot,'/models/object/',ModelFolder,'meshlab/');

FileLocationMeshItL = strcat(FileLocationWD,ModelFolder,'LBO/mesh/');
FileLocationNeighbors = strcat(FileLocationWD,ModelFolder,'neighbors/');

FileNameLBM = strcat(Model,'_LBM','_mesh_', options.dtype(1:3),'_rho',...
              num2str(options.rho),'_', options.htype,'.mat'); 
          


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[PointCloud.Location, PointCloud.Face, PointCloud.Normal, PointCloud.Signal]...
                = read_ply_all_elements( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
PointCloud = findMeshResolution(PointCloud, 'Model');



FileNameNeighbors = strcat(Model,'_Neighbors.mat');
if ~exist( strcat( FileLocationNeighbors, FileNameNeighbors), 'file')
    [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);
    save(strcat( FileLocationNeighbors, FileNameNeighbors) ,'Neighbors', '-v7.3')
else
    load(strcat( FileLocationNeighbors, FileNameNeighbors), 'Neighbors');
end


PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l_tau = set_ltau(PointCloud.Resolution);
t_0 = set_t0(PointCloud.Resolution);
NumSteps = set_NumSteps(PointCloud.Resolution, t_0);

if strcmp(options.htype, 'psp')
    options.hs = set_hs(PointCloud.Resolution);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        



% [PK1, PK2, PD1, PD2, MK, GK] = findPointCurvatures(PointCloud, NormalRotations, Neighbors.Connect);
% clear PK1 PK2 PD1 PD2 GK NeighborFaces
Quants = quantile(PointCloud.Signal, [0.25,0.5,0.75]);
MKQuant = PointCloud.Signal;
OutOfBounds = (MKQuant > (Quants(3) + 1.5*(Quants(3)-Quants(1)))) | (MKQuant < (Quants(1) - 1.5*(Quants(3)-Quants(1))));
MKQuant(OutOfBounds) = [];
stdMK = std(MKQuant);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist(strcat(FileLocationMeshItL, FileNameLBM), 'file')
    [LBM] = makeMeshLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), options);
    save(strcat(FileLocationMeshItL, FileNameLBM), 'LBM','-v7.3');
else
    load(strcat(FileLocationMeshItL, FileNameLBM), 'LBM');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[tn, tau_n, sigma_n] = findScaleStep(k, t_0, alpha, NumSteps);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusion of Mean Curvature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : length(NoiseVec)
	
    for i = 1 : NumIter
        
        sprintf('Std %0.1f : %d',NoiseVec(j),i)
        
        SignalNoisy = PointCloud.Signal + NoiseVec(j)*stdMK*randn(PointCloud.LocationCount,1);
   
        [Signal, IterCount] = performBDFDiffusion(SignalNoisy, LBM, alpha, tau_n, NumSteps, l_tau);                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find Difference of Gaussian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        DoG = buildDoG(Signal, sigma_n, DoGNormalize);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Detect Extrema
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Distance);
        Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
        Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
        Keypoint.Scale = sigma_n(Keypoint.Level);
        Keypoint.Count = length(Keypoint.Scale);
        
 
 
        NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, l_scale, l_range1, 'sigma', DoGNormalize, CompareMethod);               
        
        % SubKeypoint = findSubKeypoint(Keypoint, ScaleParameterAbsolute, DoG, PointCloud, Neighbors, NeighborFaces);
        
        
        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'SignalNoise/',Destination,'/Std_',num2str(NoiseVec(j)),'/');
        FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
        
        
        FileName = strcat('NMSKeypoint','_sigma_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
        
        
        NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, l_scale, l_range2, 'ebar', DoGNormalize, CompareMethod);
        FileName = strcat('NMSKeypoint','_ebar_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
    end
end




for i = 1
    
    sprintf('No Noise : %d',i)
    
    [Signal, IterCount] = performBDFDiffusion(PointCloud.Signal, LBM, alpha, tau_n, NumSteps, l_tau);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find Difference of Gaussian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DoG = buildDoG(Signal, sigma_n, DoGNormalize);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detect Extrema
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Distance);
    Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
    Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
    Keypoint.Scale = sigma_n(Keypoint.Level);
    Keypoint.Count = length(Keypoint.Scale);
    
    NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, l_scale, l_range1, 'sigma', DoGNormalize, CompareMethod);
 
    %     SubKeypoint = findSubKeypoint(Keypoint, ScaleParameterAbsolute, DoG, PointCloud, Neighbors.Connect, NeighborFaces.Connect);
    
    
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'SignalNoise/',Destination,'/');
    FileName = strcat('Keypoint','.mat');
    
    save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
    
    FileName = strcat('NMSKeypoint_sigma','.mat');
    save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')

    
    NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, l_scale, l_range2, 'ebar', DoGNormalize, CompareMethod);
    FileName = strcat('NMSKeypoint_ebar','.mat');
    save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
  
end








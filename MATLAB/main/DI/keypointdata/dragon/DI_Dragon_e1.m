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

porder = 4; % order of interpolation
dim = 3; % dimension
Lorder = 2; % Cartesian Laplace order
spacing = 0.01; % spacing of embedding grid


ShowPlot = 0;
Model = 'dragon/Dragon_e1_100000';

BDF = 4;
tauFraction = 1/8;
maxTauNumer = 2.5;
DoGNormalize = 1;

NumIter = 50;
NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];

t_scale = 0.7;
t_DoG = 0.9;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
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
PointCloud = findMeshNormals(PointCloud);
NormalRotations = findNormalsRotation(PointCloud.Normal);



% % % % % % % % % %
tau = spacing * tauFraction;
MaxTau = maxTauNumer / spacing;
NumSteps = round(MaxTau);
% % % % % % % % % %


[Neighbors, NeighborFaces] = findAdjacentNeighbors(PointCloud);

[PK1, PK2, PD1, PD2, MK, GK] = findPointCurvatures(PointCloud, NormalRotations, Neighbors.Connect);

PointCloud.Signal = MK;
clear PK1 PK2 PD1 PD2 GK NeighborFaces
stdMK = std(MK);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                                      
[ItL, Eplot, CP, CPFACE] = makeImplicitLaplaceBeltrami( PointCloud, porder, Lorder, spacing, dim, BDF, tau, alpha);

FaceInterpolateWeights = interpBarycenterTriangle(PointCloud, CP, CPFACE);

% CPSignal = FaceInterpolateWeights * PointCloud.Signal;

ScaleParameter = findScaleParamter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect Extrema
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j = 1 : length(NoiseVec)
    
    for i = 1 : NumIter
        i
        PointCloud.Signal = MK + NoiseVec(j)*stdMK*rand(PointCloud.LocationCount,1);
        
        CPSignal = FaceInterpolateWeights * PointCloud.Signal;
        
        CPSignalDiffused = performBDFDiffusion(CPSignal, NumSteps, ItL);

        Signal = performEplotProjection(CPSignalDiffused, Eplot);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find Difference of Gaussian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Detect Extrema
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Keypoint = findKeypoint(DoG, ScaleParameter, Neighbors.Distance, 'Old');


        NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_DoG);
        
        % SubKeypoint = findSubKeypoint(Keypoint, ScaleParameterAbsolute, DoG, PointCloud, Neighbors, NeighborFaces);
        
        
        FileLocation = strcat(ProjectRoot,'/main/DI/keypointdata/dragon/Std_',num2str(NoiseVec(j)),'/');
        FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
        %     FileName = strcat('Keypoint','.mat');
        
        save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
        
        
        FileName = strcat('NMSKeypoint','_Iter',num2str(i),'.mat');
        %     FileName = strcat('NMSKeypoint','.mat');
        
        save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
        
        
    end
    
end


























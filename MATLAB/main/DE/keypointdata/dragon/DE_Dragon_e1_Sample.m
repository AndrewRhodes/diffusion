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

options.rho = 6;
options.dtype = 'geodesic';
options.htype = 'psp';
% options.dtype = 'euclidean';

ModelFolder = 'dragon/';
Model = 'Dragon_e1';

BDF = 2;
% tauFraction = 1/10;
NumIter = 20;
NumSteps = 2000;
DoGNormalize = 'NLoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

t_scale = 0.7;
% t_DoG = 0.9;
t_range = 1/2;

SamplePercVec = [0.06, 0.055, 0.05, 0.045, 0.04];
SamplePercTrue = 0.065;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model File Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');
FileNameModelOff = strcat(ModelFolder,Model,'.off');

TmpLocation = strcat(ProjectRoot,'/models/object/',ModelFolder,'meshlab/');
if ~exist(TmpLocation,'dir')
    mkdir(TmpLocation)
end
SaveLocation = strcat(ProjectRoot, '/main/DE/keypointdata/',ModelFolder);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PointCloudOriginal.Location, PointCloudOriginal.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloudOriginal.LocationCount = size(PointCloudOriginal.Location,1);
PointCloudOriginal.FaceCount = size(PointCloudOriginal.Face,1);


for i = 1 : length(SamplePercVec)
    i
    
    % % % % Reduce the point cloud size % % % %
    PC.faces = PointCloudOriginal.Face;
    PC.vertices = PointCloudOriginal.Location;
    
    DownSampleFacesNum = round(SamplePercVec(i) * PointCloudOriginal.FaceCount);
    PC = reducepatch(PC, DownSampleFacesNum);
    TmpName = strcat(TmpLocation,Model,'_i',num2str(i),'.ply');
    ply_write(TmpName, PC.faces, PC.vertices)
    clear PC
    
    
    % % % % Run python script % % % %
    % Script open meshlab, curvature estimation, model cleaning, and saves as ply.
    
    [status,cmdout] = system("python3 "...
                    + strcat(ProjectRoot,'/src/python/callMeshlab_CleanMesh.py')...
                    + " " + TmpName + " "...
                    + strcat(Model,'_i',num2str(i),'_clean.ply'));
    
    % % % % Read .ply from meshlab % % % %            
    [PointCloud.Location, PointCloud.Face, PointCloud.Normal, PointCloud.Signal]...
        = read_ply_all_elements( strcat(TmpLocation,Model,'_i',num2str(i),'_clean.ply') );
    
    % % % % Remove clean .ply from folder % % % %            
    [status,cmdout] = system("rm " + strcat(TmpLocation,Model,'_i',num2str(i),'_clean.ply'));
    
    
    PointCloud.LocationCount = size(PointCloud.Location,1);
    PointCloud.FaceCount = size(PointCloud.Face, 1);
%     PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
    PointCloud = findMeshResolution(PointCloud, 'Model');
%     PointCloud = findMeshNormals(PointCloud);


    % % % % Save point cloud as .off and .ply % % % %
    FinalName = strcat(FileLocationModel, ModelFolder, Model, '_',num2str(PointCloud.FaceCount));
    save_off(PointCloud.Location, PointCloud.Face, strcat(FinalName, '.off') );
    ply_write(strcat(FinalName, '.ply'), PointCloud.Face, PointCloud.Location);
    
    
    
    % % % % Find and save neighbors % % % %
    NeighborsFileName = strcat(SaveLocation,Model,'_',num2str(PointCloud.FaceCount),'_Neighbors.mat');
    if ~exist(NeighborsFileName, 'file')
        [Neighbors, ~, PointCloud] = findAdjacentNeighbors(PointCloud);
        save(NeighborsFileName, 'Neighbors', '-v7.3')
    else
        load(NeighborsFileName, 'Neighbors')
        PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);
    end



    % % % % Build the mesh Laplace Beltrami % % % %
    options.hs = PointCloud.Resolution/2;
    tau = options.hs^2/4;
    
    ItLFileName = strcat(SaveLocation,Model,'_',num2str(PointCloud.FaceCount),'_ItL_rho',num2str(options.rho),'_dtype_',options.dtype,'.mat');
    
    if ~exist(ItLFileName, 'file')
        ItL = makeExplicitLaplaceBeltrami( strcat(FinalName, '.off') , options, BDF, tau, 1);
        save(ItLFileName, 'ItL', '-v7.3')
    else
        load(ItLFileName, 'ItL')
    end
    
    % % % % Find Keypoints % % % %
    ScaleParameter = findScaleParamter(tau, alpha, NumSteps, 'Laplacian', 'Natural');
    
    ParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);
    
    Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
    
    DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);
    
    Keypoint = findKeypoint(DoG, PointCloud, ScaleParameter, Neighbors.Distance, KeypointMethod, CompareMethod);

    NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);
  
    
    
    
    % % % % Save keypoints % % % %
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'Downsample/Perc_',num2str(SamplePercVec(i)),'/');
    if ~exist(FileLocation ,'dir')
        mkdir(FileLocation)
    end
%     
%    
%     FileName = strcat('Keypoint.mat');
%     save(fullfile(strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'Downsample/'), FileName), 'Keypoint', '-v7.3')
%     
%     FileName = strcat('NMSKeypoint.mat');
%     save(fullfile(strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'Downsample/'), FileName), 'NMSKeypoint', '-v7.3')
% 


    FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
    save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
    
    FileName = strcat('NMSKeypoint','_Iter',num2str(i),'.mat');
    save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')

end



















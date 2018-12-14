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

Destination = 'Umbrella';
ModelFolder = 'bunny/';
Model = 'Bunny_e1';

BDF = 1;
NumSteps = 2000;
tauFraction = 8;
DoGNormalize = 'NLoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

t_scale = 1/sqrt(2);
t_range = 1/2;

SamplePercVec = [1,0.98,0.95,0.93,0.9,0.88,0.85,0.83];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model File Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');
FileNameModelOff = strcat(ModelFolder,Model,'.off');
FileLocationWD = '/media/andrew/WDRhodes/diffusiondata/';

TmpLocation = strcat(ProjectRoot,'/models/object/',ModelFolder,'meshlab/');
if ~exist(TmpLocation,'dir')
    mkdir(TmpLocation)
end
SaveLocation = strcat(ProjectRoot, '/main/DE/keypointdata/',ModelFolder);

FileLocationMeshItL = strcat(FileLocationWD,ModelFolder,'LBO/umbrella/');
FileLocationNeighbors = strcat(FileLocationWD,ModelFolder,'neighbors/');


setTau = @(e_bar) e_bar / tauFraction;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PointCloudOriginal.Location, PointCloudOriginal.Face, PointCloudOriginal.Normal, PointCloudOriginal.Signal]...
                = read_ply_all_elements( fullfile( FileLocationModel, FileNameModelPly ) );

            
PointCloudOriginal.LocationCount = size(PointCloudOriginal.Location,1);
PointCloudOriginal.FaceCount = size(PointCloudOriginal.Face, 1);
PointCloudOriginal.FaceArea = findFaceArea(PointCloudOriginal.Location,PointCloudOriginal.Face);
PointCloudOriginal = findMeshResolution(PointCloudOriginal, 'Model');


PC.faces = PointCloudOriginal.Face;
PC.vertices = PointCloudOriginal.Location;

DownSampleFacesNum = round(SamplePercVec(end) * PointCloudOriginal.FaceCount);
PC = reducepatch(PC, DownSampleFacesNum);
PointCloud.Face = PC.faces;
PointCloud.Location = PC.vertices;
PointCloud = findMeshResolution(PointCloud, 'Model');

MaxScale = findScaleParameter(setTau(PointCloud.Resolution), alpha, NumSteps, 'Laplacian', 'Natural');
MaxScale = MaxScale(end)

clear PointCloud PC DownSampleFacesNum


for i = 1 : length(SamplePercVec) 
    
    sprintf('DownSample %0.3f %% ',SamplePercVec(i))
    
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
    
    [status,cmdout] = system("python3 " ...
                    + strcat(ProjectRoot,'/src/python/callMeshlab_CleanMesh.py') ...
                    + " " + TmpName + " "...
                    + strcat(Model,'_i',num2str(i),'_clean.ply'));
    
    % % % % Read .ply from meshlab % % % %            
    [PointCloud.Location, PointCloud.Face, PointCloud.Normal, PointCloud.Signal]...
        = read_ply_all_elements( strcat(TmpLocation,Model,'_i',num2str(i),'_clean.ply') );    
    
    PointCloud.LocationCount = size(PointCloud.Location,1);
    PointCloud.FaceCount = size(PointCloud.Face, 1);
    PointCloud = findMeshResolution(PointCloud, 'Model');   
    
    % % % % Remove clean .ply from folder % % % %
    [status,cmdout] = system("rm " + strcat(TmpLocation,Model,'_i',num2str(i),'_clean.ply'));
    
    
    % % % % Save point cloud as .off and .ply % % % %
    FileNameObject = strcat(FileLocationModel,ModelFolder,'Sample/',Model,'_',num2str(PointCloud.FaceCount));
    FileNamePly = strcat(FileNameObject, '.ply');
    FileNameOff = strcat(FileNameObject, '.off');
    
    save_off(PointCloud.Location, PointCloud.Face, FileNameOff );
    write_ply_all_elements(FileNamePly ,PointCloud.Face,PointCloud.Location,PointCloud.Normal,PointCloud.Signal)

 
    % % % % Find and save neighbors % % % %
    NeighborsFileName = strcat(FileLocationNeighbors,Model,'_',num2str(PointCloud.FaceCount),'_Neighbors.mat');
    if ~exist(NeighborsFileName, 'file')
        [Neighbors, ~, PointCloud] = findAdjacentNeighbors(PointCloud);
        save(NeighborsFileName, 'Neighbors', '-v7.3')
    else
        load(NeighborsFileName, 'Neighbors')
        PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tau = setTau(PointCloud.Resolution);
    NumSteps = round( (( MaxScale^2 ) / (2 * alpha * tau)) ) + 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ItL = makeUmbrellaLaplaceBeltrami(PointCloud, Neighbors.Connect, tau, alpha, BDF);
        
    % % % % Find Keypoints % % % %
    ScaleParameter = findScaleParameter(tau, alpha, NumSteps, 'Laplacian', 'Natural');
    
%     ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);
    
    Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
    
    DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);
    
%     Keypoint = findKeypoint(DoG, PointCloud, ScaleParameter, Neighbors.Distance, KeypointMethod, CompareMethod);
    [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Connect);
    Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
    Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
    Keypoint.Scale = ScaleParameter(Keypoint.Level);        
        
    NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);
    

    
    FileNameSubfix = strcat('Face_',num2str(num2str(PointCloud.FaceCount)),...
                      '_N',num2str(NumSteps),'.mat');
    FileNameKeypoint = strcat('Keypoint_', FileNameSubfix);
    FileNameNMSKeypoint = strcat('NMSKeypoint_', FileNameSubfix);

    if i == 1
        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'Sample/',Destination,'/');
        save(fullfile(FileLocation, FileNameKeypoint), 'Keypoint', '-v7.3')
        save(fullfile(FileLocation, FileNameNMSKeypoint), 'NMSKeypoint', '-v7.3')
    else
        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'Sample/',Destination,'/Cases/');
        save(fullfile(FileLocation, FileNameKeypoint), 'Keypoint', '-v7.3')
        save(fullfile(FileLocation, FileNameNMSKeypoint), 'NMSKeypoint', '-v7.3')
    end
    
    clear Keypoint NMSKeypoint Signal DoG ItL
end



















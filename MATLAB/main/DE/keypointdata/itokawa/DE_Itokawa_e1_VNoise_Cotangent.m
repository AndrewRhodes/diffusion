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

Destination = 'Cot';
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


FileLocationMeshItL = strcat(FileLocationWD,ModelFolder,'LBO/cotangent/');
FileLocationNeighbors = strcat(FileLocationWD,ModelFolder,'neighbors/');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PointCloudOriginal.Location, PointCloudOriginal.Face, PointCloudOriginal.Normal, PointCloudOriginal.Signal]...
                = read_ply_all_elements( fullfile( FileLocationModel, FileNameModelPly ) );


PointCloudOriginal.LocationCount = size(PointCloudOriginal.Location,1);
PointCloudOriginal.FaceCount = size(PointCloudOriginal.Face, 1);
PointCloudOriginal.FaceArea = findFaceArea(PointCloudOriginal.Location,PointCloudOriginal.Face);
PointCloudOriginal = findMeshResolution(PointCloudOriginal, 'Model');

NumSteps = set_NumSteps(PointCloudOriginal.Resolution, set_t0(PointCloudOriginal.Resolution));


for j = 1 : length(NoiseVec)

    for i = 1 : NumIter
        
        sprintf('Std %0.1f : %d',NoiseVec(j),i)
        
        FileNameObject = strcat(FileLocationModel,ModelFolder,'VertexNoise/',Model,'_sigma',num2str(NoiseVec(j)),'_iter',num2str(i));
        FileNamePly = strcat(FileNameObject, '.ply');
        FileNameOff = strcat(FileNameObject, '.off');
        
        if ~exist(FileNamePly, 'file')
        
            Location = PointCloudOriginal.Location + NoiseVec(j) * PointCloudOriginal.Resolution * ( randn(PointCloudOriginal.LocationCount,1) .* PointCloudOriginal.Normal );
            TmpName = strcat(TmpLocation,Model,'_sigma',num2str(NoiseVec(j)),'.ply');
            ply_write(TmpName, PointCloudOriginal.Face, Location);

            % % % % Run python script % % % %
            % Script open meshlab, curvature estimation,

            [status,cmdout] = system("python3 "...
                        + strcat(ProjectRoot,'/src/python/callMeshlab_findCurvature.py')...
                        + " " + TmpName + " "...
                        + strcat(Model,'_sigma',num2str(NoiseVec(j)),'_iter',num2str(i),'.ply'));

            % % % % Read .ply from meshlab % % % %            
            [PointCloud.Location, PointCloud.Face, PointCloud.Normal, PointCloud.Signal]...
                = read_ply_all_elements( FileNamePly );

            PointCloud.LocationCount = size(PointCloud.Location,1);
            PointCloud.FaceCount = size(PointCloud.Face, 1);
            PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);        
            PointCloud = findMeshResolution(PointCloud, 'Model');

            % % % % Save point cloud as .off and .ply % % % %
            save_off(PointCloud.Location, PointCloud.Face, FileNameOff );
        
        else
            [PointCloud.Location, PointCloud.Face, PointCloud.Normal, PointCloud.Signal]...
            = read_ply_all_elements( FileNamePly );
    
            PointCloud.LocationCount = size(PointCloud.Location,1);
            PointCloud.FaceCount = size(PointCloud.Face, 1);
            PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);        
            PointCloud = findMeshResolution(PointCloud, 'Model');
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup Neighbors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FileNameNeighbors = strcat(Model,'_Neighbors','_sigma',num2str(NoiseVec(j)),'_iter',num2str(i),'.mat');
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup Laplace-Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               

        FileNameLBM = strcat(Model,'_LBM','_Cot_','_sigma',num2str(NoiseVec(j)),...
                  '_iter',num2str(i),'.mat');
                                  

        if ~exist( strcat(FileLocationMeshItL, FileNameLBM),'file')
            LBM = makeCotangentLaplaceBeltrami( FileNameOff );
            save(strcat(FileLocationMeshItL, FileNameLBM),'LBM','-v7.3');
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

       [Signal, IterCount] = performBDFDiffusion(PointCloud.Signal, LBM, alpha, tau_n, NumSteps, l_tau);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find Difference of Gaussian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        DoG = buildDoG(Signal, sigma_n, DoGNormalize);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Detect Extrema
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Connect);
        Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
        Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
        Keypoint.Scale = sigma_n(Keypoint.Level);
        Keypoint.Count = length(Keypoint.Scale);
        
               
        NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, l_scale, l_range1, 'sigma', DoGNormalize, CompareMethod);               
        
        
        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'VertexNoise/',Destination,'/Std_',num2str(NoiseVec(j)),'/');
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
    
    
    FileNameNeighbors = strcat(Model,'_Neighbors.mat');
    if ~exist( strcat( FileLocationNeighbors, FileNameNeighbors), 'file')
        [Neighbors, NeighborFaces, PointCloudOriginal] = findAdjacentNeighbors(PointCloudOriginal);
        save(strcat( FileLocationNeighbors, FileNameNeighbors) ,'Neighbors', '-v7.3')
    else
        load(strcat( FileLocationNeighbors, FileNameNeighbors), 'Neighbors');
    end
    
    PointCloudOriginal = findLocalResolution(PointCloudOriginal, Neighbors.Connect);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    l_tau = set_ltau(PointCloudOriginal.Resolution);
    t_0 = set_t0(PointCloudOriginal.Resolution);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    FileNameLBM = strcat(Model,'_LBM','_Cot_','.mat'); 
          
    if ~exist(strcat(FileLocationMeshItL, FileNameLBM),'file')
        LBM = makeCotangentLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ) );         
        save(strcat(FileLocationMeshItL, FileNameLBM),'LBM','-v7.3');
    else
        load(strcat(FileLocationMeshItL, FileNameLBM), 'LBM');
    end
    
    [tn, tau_n, sigma_n] = findScaleStep(k, t_0, alpha, NumSteps);

    [Signal, IterCount] = performBDFDiffusion(PointCloudOriginal.Signal, LBM, alpha, tau_n, NumSteps, l_tau);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find Difference of Gaussian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DoG = buildDoG(Signal, sigma_n, DoGNormalize);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detect Extrema
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Connect);
    Keypoint.Normal = PointCloudOriginal.Normal(Keypoint.LocationIndex,:);
    Keypoint.Location = PointCloudOriginal.Location(Keypoint.LocationIndex,:);
    Keypoint.Scale = sigma_n(Keypoint.Level);
    Keypoint.Count = length(Keypoint.Scale);
    
    
    NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, l_scale, l_range1, 'sigma', DoGNormalize, CompareMethod);               
 
       
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'VertexNoise/',Destination,'/');
    FileName = strcat('Keypoint','.mat');
    
    save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
    
    FileName = strcat('NMSKeypoint_sigma','.mat');
    save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
    
    
    NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, l_scale, l_range2, 'ebar', DoGNormalize, CompareMethod);
    FileName = strcat('NMSKeypoint_ebar','.mat');
    save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
    
    


end





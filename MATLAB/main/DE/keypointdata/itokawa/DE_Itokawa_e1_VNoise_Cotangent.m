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


ModelFolder = 'itokawa/';
Model = 'Itokawa_e1_80000';

BDF = 1;
% tauFraction = 1/10;
NumIter = 10;
NumSteps = 2000;
DoGNormalize = 'NLoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

t_scale = 0.7;
% t_DoG = 0.9;
t_range = 1/2;
% NMSPercent = 0.1; % Percentage <= 1

NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model File Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');
FileNameModelOff = strcat(ModelFolder,Model,'.off');

TmpLocation = strcat(ProjectRoot,'/models/object/',ModelFolder,'meshlab/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PointCloudOriginal.Location, PointCloudOriginal.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );


PointCloudOriginal.LocationCount = size(PointCloudOriginal.Location,1);
PointCloudOriginal.FaceCount = size(PointCloudOriginal.Face, 1);
PointCloudOriginal.FaceArea = findFaceArea(PointCloudOriginal.Location,PointCloudOriginal.Face);
PointCloudOriginal = findMeshResolution(PointCloudOriginal, 'Model');
PointCloudOriginal = findMeshNormals(PointCloudOriginal);

load(strcat(Model,'_Neighbors.mat'),'Neighbors')
PointCloudOriginal = findLocalResolution(PointCloudOriginal, Neighbors.Connect);


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
        


        % [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);
        % save Armadillo_e1_100000_Neighbors Neighbors

        
        PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);
        % % % % % % % % % %
        tau = (PointCloudOriginal.Resolution/2)^2/4;
        % % % % % % % % % %


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup Laplace-Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FileNameItL = strcat(Model,'_ItL','_Cot','_BDF',num2str(BDF),'_sigma',num2str(NoiseVec(j)), ...
            '_iter',num2str(i),'.mat');
        
        
        if ~exist(FileNameItL,'file')
            ItL = makeCotangentLaplaceBeltrami( FileNameOff, BDF, tau, alpha);
            save(FileNameItL,'ItL','-v7.3');
        else
            load(FileNameItL, 'ItL');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Scale Parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ScaleParameter = findScaleParameter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

%         ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Diffusion of Mean Curvature
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
        Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
        
%         Signal = performBDFDiffusion_cpp(ItL, PointCloud.Signal, NumSteps);
 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find Difference of Gaussian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Detect Extrema
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %     Keypoint = findKeypoint(DoG, PointCloud, ScaleParameter, Neighbors.Distance, KeypointMethod, CompareMethod);
        [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Connect);
        Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
        Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
        Keypoint.Scale = ScaleParameter(Keypoint.Level);



        NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);


        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'VertexNoise/Cotangent/Std_',num2str(NoiseVec(j)),'/');
        FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')

        FileName = strcat('NMSKeypoint','_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
    
    
    end
end




for i = 1
    
    sprintf('No Noise : %d',i)
    
    % % % % % % % % % %
    tau = (PointCloudOriginal.Resolution/2)^2/4;
    % % % % % % % % % %
    
    FileNameItL = strcat(Model,'_ItL','_Cot','_BDF',num2str(BDF),'.mat');

    if ~exist(FileNameItL,'file')
        ItL = makeCotangentLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), BDF, tau, alpha);
        save(FileNameItL,'ItL','-v7.3');
    else
        load(FileNameItL, 'ItL');
    end
    
    
    
            
    
    load(strcat(Model,'_Curvature.mat'),'Curvature')
    MK = Curvature;
      

    ScaleParameter = findScaleParameter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

%     ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloudOriginal.Resolution);
    
    
    PointCloudOriginal.Signal = MK;
    
    Signal = performBDFDiffusion(PointCloudOriginal.Signal, NumSteps, ItL);
%     Signal = performBDFDiffusion_cpp(ItL, PointCloud.Signal, NumSteps);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find Difference of Gaussian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detect Extrema
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     Keypoint = findKeypoint(DoG, PointCloud, ScaleParameter, Neighbors.Distance, KeypointMethod, CompareMethod);
    [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Connect);
    Keypoint.Normal = PointCloudOriginal.Normal(Keypoint.LocationIndex,:);
    Keypoint.Location = PointCloudOriginal.Location(Keypoint.LocationIndex,:);
    Keypoint.Scale = ScaleParameter(Keypoint.Level);

    
    NMSKeypoint = applyNMS(PointCloudOriginal, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);

    
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'VertexNoise/Cotangent/');
    FileName = strcat('Keypoint','.mat');
    
    save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
    
    FileName = strcat('NMSKeypoint','.mat');
    save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')

end




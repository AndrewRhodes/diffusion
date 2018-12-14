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

Destination = 'Cotangent_psp_opt';
ModelFolder = 'armadillo/';
Model = 'Armadillo_e1_100000';

BDF = 1;
NumIter = 20;
NumSteps = 1000;
DoGNormalize = 'NLoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

t_scale = 1/sqrt(2);
t_range = 1/2;

NoiseVec = [0.3, 0.4, 0.5];


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


setTau = @(e_bar) e_bar^(2/5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PointCloudOriginal.Location, PointCloudOriginal.Face, PointCloudOriginal.Normal, PointCloudOriginal.Signal]...
                = read_ply_all_elements( fullfile( FileLocationModel, FileNameModelPly ) );


PointCloudOriginal.LocationCount = size(PointCloudOriginal.Location,1);
PointCloudOriginal.FaceCount = size(PointCloudOriginal.Face, 1);
PointCloudOriginal.FaceArea = findFaceArea(PointCloudOriginal.Location,PointCloudOriginal.Face);
PointCloudOriginal = findMeshResolution(PointCloudOriginal, 'Model');



for j = 1 : length(NoiseVec)
    if j == 1 
        List = 7 : NumIter;
    else
        List = 1 : NumIter;
    end
    for i = List
        
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
        tau = setTau(PointCloud.Resolution);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup Laplace-Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FileNameItL = strcat(Model,'_ItL','_Cot','_BDF',num2str(BDF),...
            '_te0.4','_a',num2str(alpha),'_sigma',num2str(NoiseVec(j)),'_iter',num2str(i),'.mat');
                
        if ~exist( strcat(FileLocationMeshItL, FileNameItL),'file')
            ItL = makeCotangentLaplaceBeltrami( FileNameOff, BDF, tau, alpha);
            save(strcat(FileLocationMeshItL, FileNameItL),'ItL','-v7.3');
        else
            load(strcat(FileLocationMeshItL, FileNameItL), 'ItL');
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


        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'VertexNoise/',Destination,'/Std_',num2str(NoiseVec(j)),'/');
        FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')

        FileName = strcat('NMSKeypoint','_Iter',num2str(i),'.mat');
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
    tau = setTau(PointCloudOriginal.Resolution);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    FileNameItL = strcat(Model,'_ItL','_Cot','_BDF',num2str(BDF),...
            '_te0.4','_a',num2str(alpha),'.mat');

    if ~exist(strcat(FileLocationMeshItL, FileNameItL),'file')
        ItL = makeCotangentLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), BDF, tau, alpha);
        save(strcat(FileLocationMeshItL, FileNameItL),'ItL','-v7.3');
    else
        load(strcat(FileLocationMeshItL, FileNameItL), 'ItL');
    end
    

    ScaleParameter = findScaleParameter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

%     ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloudOriginal.Resolution);
    
        
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

    
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'VertexNoise/',Destination,'/');
    FileName = strcat('Keypoint','.mat');
    
    save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
    
    FileName = strcat('NMSKeypoint','.mat');
    save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call Other Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DE_Armadillo_e1

DE_Armadillo_e1_Umbrella

DE_Armadillo_e1_Cotangent




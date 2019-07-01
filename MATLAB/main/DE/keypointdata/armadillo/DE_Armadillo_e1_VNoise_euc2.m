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

Destination = 'Mesh_rho4_ddr_euc'
ModelFolder = 'armadillo/';
Model = 'Armadillo_e1_100000';


BDF = 1;
NumIter = 20;
NumSteps = 2000;
DoGNormalize = 'NLoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

t_scale = 1/sqrt(2);
t_range1 = 1/2;
t_range2 = 2;

NoiseVec = [0.4];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model File Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileLocationWD = '/media/andrew/WDRhodes/diffusiondata/';
FileNameModelPly = strcat(ModelFolder,Model,'.ply');
FileNameModelOff = strcat(ModelFolder,Model,'.off');

TmpLocation = strcat(ProjectRoot,'/models/object/',ModelFolder,'meshlab/');

FileLocationMeshItL = strcat(FileLocationWD,ModelFolder,'LBO/mesh/');
FileLocationNeighbors = strcat(FileLocationWD,ModelFolder,'neighbors/');

FileNameItL = strcat(Model,'_ItL','_BDF',num2str(BDF),'_rho',...
              num2str(options.rho),'_',options.dtype(1:3),'_',...
              options.htype,'_t0.25e','_a',num2str(alpha),'.mat');
          
          
setTau = @(e_bar) 0.25*e_bar;
setHs = @(e_bar) 2*e_bar^(1/5);

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

    for i = 19 : 20
        
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
        if strcmp(options.htype, 'psp')
            options.hs = setHs(PointCloud.Resolution);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup Laplace-Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        FileNameItL = strcat(Model,'_ItL','_BDF',num2str(BDF),'_rho',...
              num2str(options.rho),'_',options.dtype(1:3),'_',...
              options.htype,'_t0.25e','_a',num2str(alpha),...
              '_sigma',num2str(NoiseVec(j)), '_iter',num2str(i),'.mat');
          
        
        if ~exist( strcat(FileLocationMeshItL, FileNameItL),'file')
            ItL = makeMeshLaplaceBeltrami( FileNameOff, options, BDF, tau, alpha);
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
        [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Distance);
        Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
        Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
        Keypoint.Scale = ScaleParameter(Keypoint.Level);

        NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range1, 'sigma', DoGNormalize, CompareMethod);     


        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'VertexNoise/',Destination,'/Std_',num2str(NoiseVec(j)),'/');
        FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')

        FileName = strcat('NMSKeypoint','_sigma_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')      
    
        
        NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range2, 'ebar', DoGNormalize, CompareMethod);
        FileName = strcat('NMSKeypoint','_ebar_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
        
        
    end
end







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


alpha = 0.5;


ModelFolder = 'armadillo/';
Model = 'Armadillo_e1_100000';

BDF = 1;
% tauFraction = 1/10;
NumIter = 20;
NumSteps = 2000;
DoGNormalize = 'NLoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

t_scale = 0.7;
t_DoG = 0.9;
t_range = 1/2;

NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];

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


PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
PointCloud = findMeshResolution(PointCloud, 'Model');
PointCloud = findMeshNormals(PointCloud);


% load('Neighbors.mat')
% [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);
% save Armadillo_e1_100000_Neighbors Neighbors

load(strcat(Model,'_Neighbors.mat'),'Neighbors')
PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);


% % % % % % % % % %
% options.hs = PointCloud.Resolution/2;
% tau = options.hs^2/4;
tau = (PointCloud.Resolution/2)^2/4;
% % % % % % % % % %

load(strcat(Model,'_Curvature.mat'),'Curvature')
MK = Curvature;



% [PK1, PK2, PD1, PD2, MK, GK] = findPointCurvatures(PointCloud, NormalRotations, Neighbors.Connect);
% clear PK1 PK2 PD1 PD2 GK NeighborFaces
Quants = quantile(MK, [0.25,0.5,0.75]);
MKQuant = MK;
OutOfBounds = (MKQuant > (Quants(3) + 1.5*(Quants(3)-Quants(1)))) | (MKQuant < (Quants(1) - 1.5*(Quants(3)-Quants(1))));
MKQuant(OutOfBounds) = [];
stdMK = std(MKQuant);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ItL = makeUmbrellaLaplaceBeltrami(PointCloud, Neighbors.Connect, tau, alpha, BDF);
% save( strcat(Model,'_ItL_umb_BDF',num2str(BDF),'.mat') , '-v7.3')
% load( strcat(Model,'_ItL_umb_BDF',num2str(BDF),'.mat'), 'ItL')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ScaleParameter = findScaleParameter(tau, alpha, NumSteps, 'Laplacian', 'Natural');
% ScaleParameter2 = findScaleParameter(sqrt(tau), alpha, NumSteps, 'Gaussian', 'Natural');

ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusion of Mean Curvature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : length(NoiseVec)

    for i = 1 : NumIter
       
        i
        PointCloud.Signal = MK + NoiseVec(j)*stdMK*randn(PointCloud.LocationCount,1);
        
        
%         Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);       
 
        Signal = performBDFDiffusion_cpp(ItL, PointCloud.Signal, NumSteps);
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find Difference of Gaussian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        DoG = buildDoG(Signal(:,1:NumSteps), ScaleParameter, DoGNormalize);
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Detect Extrema
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Should I use Neighbors.Distance or Neighbors.Connect
%         Keypoint = findKeypoint(DoG, PointCloud, ScaleParameter, Neighbors.Connect, KeypointMethod, CompareMethod);

        [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Connect);
        Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
        Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
        Keypoint.Scale = ScaleParameter(Keypoint.Level);

   
        
        NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);
        
        
        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'SignalNoise/Umbrella/Std_',num2str(NoiseVec(j)),'/');
        FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
        
        
        FileName = strcat('NMSKeypoint','_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
        
        
    end
end


% sunvector = [0;0;-1];
% sunvector = sunvector / norm(sunvector);
% intensity = 0.9*PointCloud.Normal*sunvector;
% intensity(intensity<0) = 0;
% 
% 
% figure
% colormap('gray')
% trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), intensity, 'EdgeColor', 'none');
% view(180, -90)
% axis off 
% axis equal
% hold on
% 
% [SphereX, SphereY, SphereZ] = sphere(100);
% 
% SphereX = reshape(SphereX,[],1);
% SphereY = reshape(SphereY,[],1);
% SphereZ = reshape(SphereZ,[],1);
% 
% 
% for i = 1 : NMSKeypoint.Count
%     if NMSKeypoint.Scale(i) > 5
%     SphereAtPoint = bsxfun(@plus, NMSKeypoint.Scale(i)*[SphereX, SphereY, SphereZ], NMSKeypoint.Location(i, :));
% 
%     hh = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
%     set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)
%     drawnow
%     end
% %     pause
% end





for i = 1
    i
    PointCloud.Signal = MK;
    
%     Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
    
    Signal = performBDFDiffusion_cpp(ItL, PointCloud.Signal, NumSteps);    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find Difference of Gaussian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DoG = buildDoG(Signal(:,1:NumSteps), ScaleParameter, DoGNormalize);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detect Extrema
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%    Keypoint = findKeypoint(DoG, PointCloud, ScaleParameter, Neighbors.Connect, KeypointMethod, CompareMethod);
    [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Connect);
    Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
    Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
    Keypoint.Scale = ScaleParameter(Keypoint.Level);
    
    NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);

       
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'SignalNoise/Umbrella/');
    FileName = strcat('Keypoint','.mat');
    
    save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
    
    FileName = strcat('NMSKeypoint','.mat');
    save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
    
end








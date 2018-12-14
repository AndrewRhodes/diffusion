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

options.rho = 5;
options.dtype = 'geodesic';
% options.dtype = 'euclidean';

ShowPlot = 0;
ModelFolder = 'armadillo/';
Model = 'Armadillo_e1_100000';

BDF = 2;
tauFraction = 1/10;
tauNumerator = 3000;
DoGNormalize = 'DoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

t_scale = 0.7;
t_DoG = 0.9;
t_range = 3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model File Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');
FileNameModelOff = strcat(ModelFolder,Model,'.off');

load(strcat(Model,'_Curvature.mat'),'Curvature')
MK = Curvature;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );


PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
PointCloud = findMeshResolution(PointCloud, 'Model');
PointCloud = findMeshNormals(PointCloud);


% % % % % % % % % %
tau = PointCloud.Resolution * tauFraction;
MaxTau = tauNumerator / PointCloud.Resolution;
NumSteps = tauNumerator;
% % % % % % % % % %


% [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);
load(strcat(Model,'_Neighbors.mat'),'Neighbors')
PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
% ItL = makeExplicitLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), options, BDF, tau, alpha);

load(strcat(Model,'_ItL.mat'),'ItL')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ScaleParameter = findScaleParamter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusion of Mean Curvature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PointCloud.Signal = MK;


% Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
% save(strcat(Model,'_Signal_N',num2str(NumSteps),'.mat'),'Signal', '-v7.3')

load(strcat(Model,'_Signal_N',num2str(NumSteps)),'Signal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Difference of Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);
AbsDoG = buildDoG(Signal, ScaleParameter, 'AbsDoG');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect Extrema
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load(strcat(ProjectRoot,'/main/DE/keypointdata/armadillo/SignalNoise/LongRun/Std_0.1/DoGMaximum.mat'))
% load(strcat(ProjectRoot,'/main/DE/keypointdata/armadillo/SignalNoise/LongRun/Std_0.1/Keypoint.mat'))
% load(strcat(ProjectRoot,'/main/DE/keypointdata/armadillo/SignalNoise/LongRun/Std_0.1/NMSKeypoint.mat'))
% [Keypoint, DoGMaximum] = findKeypointDoGMax(DoG, ScaleParameter, PointCloud, Neighbors.Connect, KeypointMethod, CompareMethod);

[Keypoint, DoGMaximum] = findKeypointDoGMax(DoG, ScaleParameter, PointCloud, Neighbors.Connect, KeypointMethod, 'special');
        
NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);

% save KeypointSpecial Keypoint

% NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);
        

NewKeypoint = checkKeypointSphere(PointCloud, ScaleParameter, DoG, Keypoint, 1, DoGMaximum);


% save Keypoint_T3000_DoG Keypoint
% save NMSKeypoint_T3000_DoG NMSKeypoint
% save DoGMaximum_T3000_DoG DoGMaximum

% Keypoint = findKeypoint(DoG, PointCloud, ScaleParameter, Neighbors.Connect, KeypointMethod, 'special');

% load('Keypoint_T3000_DoG.mat','Keypoint')
% load('NMSKeypoint_T3000_DoG.mat','NMSKeypoint')
% load('DoGMaximum_T3000_DoG.mat','DoGMaximum')

load('DoGMaximum_T3000_AbsDoG_MH.mat','DoGMaximum')


% Remove keypoints found before \bar{sigma}_i / \bar{e}_i < 4
% ScaleLogic = NMSKeypoint.Scale ./ PointCloud.ResolutionLocal(NMSKeypoint.Location) < 4;
% NMSKeypoint.Scale(ScaleLogic) = [];
% NMSKeypoint.Level(ScaleLogic) = [];
% NMSKeypoint.Location(ScaleLogic) = [];



ScaleLogic = Keypoint.Scale ./ PointCloud.ResolutionLocal(Keypoint.Location) < 4;
% ScaleLogic = Keypoint.Level < 100;
Keypoint.Scale(ScaleLogic) = [];
Keypoint.Level(ScaleLogic) = [];
Keypoint.Location(ScaleLogic) = [];



% % % % SubKeypoint = findSubKeypoint(Keypoint, ScaleParameterAbsolute, DoG, PointCloud, Neighbors, NeighborFaces);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the image for the video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sunvector = [0;0;-1];
% sunvector = sunvector / norm(sunvector);
% intensity = 0.9*PointCloud.Normal*sunvector;
% intensity(intensity<0) = 0;


clear vidObj
vidObj = VideoWriter('Armadillo_AbsDoGMH_DoGExtrema_NewKeypoints_Spheres.avi');
vidObj.FrameRate = 20;
open(vidObj);
[SphereX, SphereY, SphereZ] = sphere(100);
% 
SphereX = reshape(SphereX,[],1);
SphereY = reshape(SphereY,[],1);
SphereZ = reshape(SphereZ,[],1);
Min = 0;
fig = figure('Units', 'Normalized', 'OuterPosition', [0 0 0.98 0.98]);
for i = 1 : tauNumerator-1

    i
   
%     figure
    trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), Signal(:,i), 'EdgeColor', 'none');
    view(180, -90)
	ax = gca;
    axis equal
    axis off
    hold on
	cbar = colorbar('FontSize',15);
   

    
%     CurrentKeypoints = NMSKeypoint.Location(NMSKeypoint.Level==i,1);
%     plot3(PointCloud.Location(CurrentKeypoints,1), PointCloud.Location(CurrentKeypoints,2),...
%         PointCloud.Location(CurrentKeypoints,3), 'k.','Markersize',30)


    CurrentDoGMax = DoGMaximum.Location(DoGMaximum.Level==i,1);
    plot3(PointCloud.Location(CurrentDoGMax,1), PointCloud.Location(CurrentDoGMax,2),...
          PointCloud.Location(CurrentDoGMax,3), 'k.','MarkerSize', 10)

%     CurrentKeypoints = Keypoint.Location(Keypoint.Level ==i,1);
%     CurrentLowKeypoints = CurrentKeypoints( Keypoint.Sign(Keypoint.Level == i) < 0 );
%     CurrentHighKeypoints = CurrentKeypoints( Keypoint.Sign(Keypoint.Level == i) > 0 );
%     
%     plot3(PointCloud.Location(CurrentLowKeypoints,1), PointCloud.Location(CurrentLowKeypoints,2),...
%           PointCloud.Location(CurrentLowKeypoints,3), 'r.','MarkerSize', 10)
%       
%     plot3(PointCloud.Location(CurrentHighKeypoints,1), PointCloud.Location(CurrentHighKeypoints,2),...
%           PointCloud.Location(CurrentHighKeypoints,3), 'b.','MarkerSize', 10)
      
      
      
    CurrentKeypoints = Keypoint.Location;%(NewKeypoint.Level ==i);
    CurrentScales = Keypoint.Scale;%(NewKeypoint.Level == i);
 
    for j = 1 : length(CurrentKeypoints)
      SphereAtPoint = bsxfun(@plus, CurrentScales(j)*[SphereX, SphereY, SphereZ], PointCloud.Location(CurrentKeypoints(j),:));
      hh = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
      set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)
    end



%         CurrentKeypoints = NMSKeypoint.Location(1:100);
%         CurrentScales = NMSKeypoint.Scale(1:100);
%       CurrentKeypoints = NMSKeypoint.Location(NMSKeypoint.Level<89,1);
%       CurrentScales = NMSKeypoint.Scale(NMSKeypoint.Level<89,1);
% 
%       plot3(PointCloud.Location(CurrentKeypoints,1), PointCloud.Location(CurrentKeypoints,2),...
%             PointCloud.Location(CurrentKeypoints,3), 'k.','Markersize',12)
%         
%       for j = 1 : length(CurrentKeypoints)
%           SphereAtPoint = bsxfun(@plus, CurrentScales(j)*[SphereX, SphereY, SphereZ], PointCloud.Location(CurrentKeypoints(j),:));
%           hh = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
%           set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)
%       end
    

    ylabel(cbar, 'Mean Curvature', 'FontSize',30);
    ax.XLim = [-240 240];
    ax.YLim = [-115 200];
    title(sprintf('Diffusion Step: %d, Absolute Scale: %0.2f',i, ScaleParameterAbsolute(i+1)),'fontsize',25)
    caxis([-0.15,0.2])
    drawnow
    pause(0.1)
    
%     FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/armadillo/Video/Images/Extrema/ArmadilloStep_',num2str(i));
%     savefig(fig, FileLocation)
        
    writeVideo(vidObj, getframe(fig));
    
    clf(fig);
%     close(fig);
    
end
close(fig);

close(vidObj)




    






















    

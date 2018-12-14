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

ShowPlot = 0;
Model = 'armadillo/Armadillo_e1_100000';

BDF = 2;
tauFraction = 1/10;
tauNumerator = 3000;
DoGNormalize = 'NLoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

t_scale = 0.7;
t_DoG = 0.9;
t_range = 3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model File Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(Model,'.ply');
FileNameModelOff = strcat(Model,'.off');

load('ArmadilloCurvature_e1_100000.mat')
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
load('Armadillo_e1_100000_Neighbors.mat')
PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save ArmadilloItL_e1_100000 ItL
 
% ItL = makeExplicitLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), options, BDF, tau, alpha);

load('ArmadilloItL_e1_100000.mat')

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
% save Signal_T3000_DoG Signal

load('Signal_T3000_DoG.mat', 'Signal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Difference of Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect Extrema
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load(strcat(ProjectRoot,'/main/DE/keypointdata/armadillo/SignalNoise/LongRun/Std_0.1/DoGMaximum.mat'))
% load(strcat(ProjectRoot,'/main/DE/keypointdata/armadillo/SignalNoise/LongRun/Std_0.1/Keypoint.mat'))
% load(strcat(ProjectRoot,'/main/DE/keypointdata/armadillo/SignalNoise/LongRun/Std_0.1/NMSKeypoint.mat'))


% [Keypoint, DoGMaximum] = findKeypointDoGMax(DoG, ScaleParameter, PointCloud, Neighbors.Distance, KeypointMethod, CompareMethod);
%         
% NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);
        

% save Keypoint_T3000_AbsDoG_MH Keypoint
% save NMSKeypoint_T000_AbsDoG_MH NMSKeypoint
% save DoGMaximum_T3000_AbsDoG_MH DoGMaximum


load('Keypoint_T3000_AbsDoG_MH.mat','Keypoint')
load('NMSKeypoint_T000_AbsDoG_MH.mat','NMSKeypoint')
load('DoGMaximum_T3000_AbsDoG_MH.mat','DoGMaximum')

t_quant = 0.1;
DoGMaximum = findMHEdge(DoG, ScaleParameter, PointCloud, Neighbors.Distance, t_quant);




% % % % SubKeypoint = findSubKeypoint(Keypoint, ScaleParameterAbsolute, DoG, PointCloud, Neighbors, NeighborFaces);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the image for the video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sunvector = [0;0;-1];
% sunvector = sunvector / norm(sunvector);
% intensity = 0.9*PointCloud.Normal*sunvector;
% intensity(intensity<0) = 0;


clear vidObj
vidObj = VideoWriter('Armadillo_DoG_MH_ZeroCrossing_Sphere.avi');
vidObj.FrameRate = 20;
open(vidObj);
[SphereX, SphereY, SphereZ] = sphere(100);
% 
SphereX = reshape(SphereX,[],1);
SphereY = reshape(SphereY,[],1);
SphereZ = reshape(SphereZ,[],1);
Min = 0;
figMH = figure('Units', 'Normalized', 'OuterPosition', [0 0 0.98 0.98]);
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
    CurrentDoGMaxSacle = DoGMaximum.Scale(DoGMaximum.Level==i,1);
%     plot3(PointCloud.Location(CurrentDoGMax,1), PointCloud.Location(CurrentDoGMax,2),...
%           PointCloud.Location(CurrentDoGMax,3), 'k.','MarkerSize', 12)


%       CurrentKeypoints = NMSKeypoint.Location(NMSKeypoint.Level==i,1);
%       CurrentScales = NMSKeypoint.Scale(NMSKeypoint.Level==i,1);

%       plot3(PointCloud.Location(CurrentKeypoints,1), PointCloud.Location(CurrentKeypoints,2),...
%             PointCloud.Location(CurrentKeypoints,3), 'k.','Markersize',50)
        
      for j = 1 : length(CurrentDoGMax)
          SphereAtPoint = bsxfun(@plus, CurrentDoGMaxSacle(j)*[SphereX, SphereY, SphereZ], PointCloud.Location(CurrentDoGMax(j),:));
          hh = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
          set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)
      end
    

    ylabel(cbar, 'Mean Curvature', 'FontSize',30);
    ax.XLim = [-240 240];
    ax.YLim = [-115 200];
    title(sprintf('Diffusion Step: %d, Absolute Scale: %0.2f',i, ScaleParameterAbsolute(i+1)),'fontsize',25)
    caxis([-0.15,0.2])
    drawnow
    pause(0.1)
    
%     FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/armadillo/Video/Images/Extrema/ArmadilloStep_',num2str(i));
%     savefig(figMH, FileLocation)
        
    writeVideo(vidObj, getframe(figMH));
    
    clf(figMH);
    
end


close(vidObj)




    






















    

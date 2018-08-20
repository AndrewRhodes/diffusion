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
Model = 'armadillo/Armadillo_e1_100000';

BDF = 2;
tauFraction = 1/10;
tauNumerator = 3000;
DoGNormalize = 'DoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

t_scale = 0.7;
t_DoG = 0.9;


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

% PointCloud.Location = PointCloud.Location.*1.9868;

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
PointCloud = findMeshResolution(PointCloud, 'Model');
PointCloud = findMeshNormals(PointCloud)
NormalRotations = findNormalsRotation(PointCloud.Normal);

% ply_write('Armadillo_e1', PointCloud.Face, PointCloud.Location)
% save_off(PointCloud.Location, PointCloud.Face, 'Armadillo_e1.off')


% % % % % % % % % %
tau = PointCloud.Resolution * tauFraction;
MaxTau = tauNumerator / PointCloud.Resolution;
NumSteps = round(MaxTau);
NumSteps = tauNumerator;
% % % % % % % % % %


% [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);
load('Armadillo_e1_100000_Neighbors.mat')
PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);



% [PK1, PK2, PD1, PD2, MK, GK] = findPointCurvatures(PointCloud, NormalRotations, Neighbors.Connect);
% clear PK1 PK2 PD1 PD2 GK NeighborFaces
stdMK = std(MK);


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


Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Difference of Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect Extrema
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Keypoint, DoGMaximum] = findKeypointDoGMax(DoG, ScaleParameter, PointCloud, Neighbors.Distance, KeypointMethod, CompareMethod);
        
NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);
        
% SubKeypoint = findSubKeypoint(Keypoint, ScaleParameterAbsolute, DoG, PointCloud, Neighbors, NeighborFaces);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the image for the video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sunvector = [0;0;-1];
% sunvector = sunvector / norm(sunvector);
% intensity = 0.9*PointCloud.Normal*sunvector;
% intensity(intensity<0) = 0;


clear vidObj
% vidObj = VideoWriter('ArmadilloVideoMeshLabMaxNLoGNMSSphere.avi');
% vidObj.FrameRate = 20;
% open(vidObj);
[SphereX, SphereY, SphereZ] = sphere(100);
% 
SphereX = reshape(SphereX,[],1);
SphereY = reshape(SphereY,[],1);
SphereZ = reshape(SphereZ,[],1);
Min = 250;
for i = 1 : 10 : tauNumerator-1

    i
    fig = figure('Units', 'Normalized', 'OuterPosition', [0 0 0.98 0.98]);
%     figure
    trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), Signal(:,i), 'EdgeColor', 'none');
    view(180, -90)
	ax = gca;
    axis equal
    axis off
    hold on
	cbar = colorbar('FontSize',15);
    caxis([-0.15,0.2])

    CurrentKeypoints = NMSKeypoint.Location(NMSKeypoint.Level==i,1);
    plot3(PointCloud.Location(CurrentKeypoints,1), PointCloud.Location(CurrentKeypoints,2),...
        PointCloud.Location(CurrentKeypoints,3), 'k.','Markersize',10)


%     CurrentDoGMax = DoGMaximum.Location(DoGMaximum.Level==i,1);
%     plot3(PointCloud.Location(CurrentDoGMax,1), PointCloud.Location(CurrentDoGMax,2),...
%           PointCloud.Location(CurrentDoGMax,3), 'k.','MarkerSize', 10)


%       CurrentKeypoints = NMSKeypoint.Location(NMSKeypoint.Level==i,1);
%       CurrentScales = NMSKeypoint.Scale(NMSKeypoint.Level==i,1);

%       plot3(PointCloud.Location(CurrentKeypoints,1), PointCloud.Location(CurrentKeypoints,2),...
%             PointCloud.Location(CurrentKeypoints,3), 'k.','Markersize',50)
        
%       for j = 1 : length(CurrentKeypoints)
%           SphereAtPoint = bsxfun(@plus, CurrentScales(j)*[SphereX, SphereY, SphereZ], PointCloud.Location(CurrentKeypoints(j),:));
%           hh = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
%           set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)
%       end
    

    ylabel(cbar, 'Mean Curvature', 'FontSize',30);
    ax.XLim = [-240 240];
    ax.YLim = [-115 200];
    title(sprintf('Diffusion Step: %d, Absolute Scale: %0.2f',i, ScaleParameterAbsolute(i+1)),'fontsize',25)
    drawnow
    pause(1)
    
%     FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/armadillo/Images/Second/ArmadilloStep_',num2str(i));
%     savefig(fig, FileLocation)
        
%     writeVideo(vidObj, getframe(fig));
    
    close(fig);
    
end


% close(vidObj)






% vidObj = VideoWriter('ArmadilloVideo.avi');
% vidObj.FrameRate = 10;
% open(vidObj);
% 
% 
% Min = 1000;
% for i = Min : tauNumerator
%     
%     FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/armadillo/Images/First/ArmadilloStep_',num2str(i));
%     
%     openfig(FileLocation)
%     ax = gca;
%     ax.XLim = [-240 240];
%     ax.YLim = [-115 200];
% 
%     writeVideo(vidObj, getframe(gca));
%     
%     close(gcf);
% end
%  
% close(vidObj)
 
 

% FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/armadillo/LongRun/Std_',num2str(NoiseVec(j)),'/');
% FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');        
% save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
% 
% 
% FileName = strcat('NMSKeypoint','_Iter',num2str(i),'.mat');
% save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')


    






















    

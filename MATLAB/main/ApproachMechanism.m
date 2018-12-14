


close all
clear
clc

global ProjectRoot; % Additional Paths
addpath(genpath('~/Desktop/MLIDAR-1.0'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


alpha = 1/2;

options.rho = 4;
options.dtype = 'geodesic'; % 'euclidean', 'geodesic'
options.htype = 'ddr'; % 'psp', 'ddr'



ModelFolder = 'itokawa/';
Model = 'Itokawa_e1_80000';


BDF = 1;
tauFraction = 8;
NumSteps = 2000;

% KeypointSearchMethod = 'Global'; % 'Local', 'Global'
% Radius = 1;

DoGNormalize = 'NLoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

t_scale = 0.7;
t_dist = 0.5;
t_gc = 0.25;
t_corr = 0.25;

% t_DoG = 0.9;
t_range = 1/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model File Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');
FileNameModelOff = strcat(ModelFolder,Model,'.off');
FileLocationWD = '/media/andrew/WDRhodes/diffusiondata/';

TmpLocation = strcat(ProjectRoot,'/models/object/',ModelFolder,'meshlab/');

FileLocationMeshItL = strcat(FileLocationWD,ModelFolder,'LBO/mesh/');
FileLocationNeighbors = strcat(FileLocationWD,ModelFolder,'neighbors/');

setTau = @(e_bar) e_bar/tauFraction;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ModelPointCloud.Location, ModelPointCloud.Face, ModelPointCloud.Normal, ModelPointCloud.Signal]...
                = read_ply_all_elements( fullfile( FileLocationModel, FileNameModelPly ) );

ModelPointCloud.LocationCount = size(ModelPointCloud.Location,1);
ModelPointCloud.FaceCount = size(ModelPointCloud.Face, 1);
ModelPointCloud.FaceArea = findFaceArea(ModelPointCloud.Location,ModelPointCloud.Face);
ModelPointCloud = findMeshResolution(ModelPointCloud, 'Model');
% ModelPointCloud = findMeshNormals(ModelPointCloud);


% [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);


FileNameNeighbors = strcat(Model,'_Neighbors.mat');
if ~exist( strcat( FileLocationNeighbors, FileNameNeighbors), 'file')
    [ModelNeighbors, ModelNeighborFaces, ModelPointCloud] = findAdjacentNeighbors(ModelPointCloud);
    save(strcat( FileLocationNeighbors, FileNameNeighbors) ,'Neighbors', '-v7.3')
else
    load(strcat( FileLocationNeighbors, FileNameNeighbors), 'Neighbors');
    ModelNeighbors = Neighbors;
    clear Neighbors;
end

ModelPointCloud = findLocalResolution(ModelPointCloud, ModelNeighbors.Connect);
% options.hs = PointCloud.Resolution/2;

% load(strcat(Model,'_Curvature.mat'),'Curvature')
% ModelPointCloud.Signal = Curvature;

% options.hs = ModelPointCloud.Resolution/2;
% tau = options.hs^2/4;
tau = setTau(ModelPointCloud.Resolution);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load('NMSKeypoint.mat')
% ModelKeypoint = NMSKeypoint;
% clear NMSKeypoint


FileNameItL = strcat(Model,'_ItL','_BDF',num2str(BDF),'_rho',...
                    num2str(options.rho),'_',options.dtype,'_',...
                    options.htype,'_ebar-',num2str(tauFraction),'.mat');
                               
if ~exist(strcat(FileLocationMeshItL, FileNameItL), 'file')
    ItL = makeMeshLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), options, BDF, tau, alpha);
    save(strcat(FileLocationMeshItL, FileNameItL),'ItL','-v7.3');
else
    load(strcat(FileLocationMeshItL, FileNameItL), 'ItL');
end

ModelScaleParameter = findScaleParameter(tau, alpha, NumSteps, 'Laplacian', 'Natural');
% ModelScaleParameterAbsolute = bsxfun(@plus, ModelScaleParameter, ModelPointCloud.Resolution);

ModelSignal = performBDFDiffusion(ModelPointCloud.Signal, NumSteps, ItL);

ModelDoG = buildDoG(ModelSignal, ModelScaleParameter, DoGNormalize);

[ModelKeypoint.LocationIndex, ModelKeypoint.Level] = findKeypoint_c(ModelDoG, ModelNeighbors.Distance);
ModelKeypoint.Normal = ModelPointCloud.Normal(ModelKeypoint.LocationIndex,:);
ModelKeypoint.Location = ModelPointCloud.Location(ModelKeypoint.LocationIndex,:);
ModelKeypoint.Scale = ModelScaleParameter(ModelKeypoint.Level);
        
% ModelKeypoint = findKeypoint(ModelDoG, ModelPointCloud, ModelScaleParameter, ModelNeighbors.Distance, KeypointMethod, CompareMethod);
% ModelKeypoint.Scale = ModelKeypoint.Scale + ModelPointCloud.ResolutionLocal(ModelKeypoint.LocationIndex);


ModelNMSKeypoint = applyNMS(ModelPointCloud, ModelDoG, ModelKeypoint, t_scale, t_range, DoGNormalize, CompareMethod);


ModelFeatures = buildHistogram(ModelPointCloud, ModelNMSKeypoint);

% load('Features.mat','Features')
% ModelFeatures = Features;
% clear Features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take 3D image of model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LOSVectors = importdata('/media/andrew/WDRhodes/GitProjects/pose/MATLAB_PointCloudDescriptors/OURCVFH/LOS/LOSVectors256_20.mat');


[~, octreeRes] = centerPoint(ModelPointCloud.Location, ModelPointCloud.Face);
AdjList = generate_adjacent_face_indices(ModelPointCloud.Face);


v = [0,2,1];
v = v/norm(v);
e = 25*pi/180;
Rot = v_to_R(e*v);
Rot = eye(3,3);
trans = [0,0,1000];

verticesInCamera = bsxfun(@plus, (Rot * ModelPointCloud.Location')', trans);

[PointDistances,~] = mex_pcs(LOSVectors, verticesInCamera(:,1:3), ModelPointCloud.Face, AdjList, octreeRes*3);

Points = bsxfun(@times, LOSVectors(PointDistances>0,:), PointDistances(PointDistances>0));
ScenePointCloud.Location = Points;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Point Cloud LBO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt.rho = 4;
ScenePointCloud.Resolution = findMeshResolution(ScenePointCloud, 'Scene');
PCD_LBO = pcdlp_matrix(Points, 3, opt);
tauPC = settau(ScenePointCloud.Resolution);

PCItL = eye(size(PCD_LBO)) - alpha*tauPC*PCD_LBO;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to find edge 'pixels' around object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Then remove keypoints if they are closer to the edge is within their
% scale support size.







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Save 3D Points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saves as a .ply file with only vertices. No faces.

if ~exist(TmpLocation,'dir')
    mkdir(TmpLocation)
end

TmpName = strcat(TmpLocation,'Points3D.ply');
write_pointcloud_ply(TmpName, Points);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Run python script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script open meshlab, applies normal construction, mesh reconstruction,
% curvature estimation, model cleaning, and saves as ply.


[status,cmdout] = system("python3 " + strcat(ProjectRoot,'/src/python/callMeshlab.py') + " " + TmpName + " " + 'Mesh3D.ply');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Load New mesh model from Meshlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[ScenePointCloud.Location, ScenePointCloud.Face, ScenePointCloud.Normal, ScenePointCloud.Signal] = read_ply_all_elements( strcat(ProjectRoot,'/models/object/itokawa/meshlab/Mesh3D.ply') );

save_off(ScenePointCloud.Location, ScenePointCloud.Face, strcat(ProjectRoot,'/models/object/itokawa/meshlab/Mesh3D.off'));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Construct Scale-Space and find keypoints, features. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ScenePointCloud = findMeshResolution(ScenePointCloud, 'Model');
ScenePointCloud.LocationCount = size(ScenePointCloud.Location,1);
ScenePointCloud.FaceCount = size(ScenePointCloud.Face,1);

[SceneNeighbors, ~, ~] = findAdjacentNeighbors(ScenePointCloud);
ScenePointCloud = findLocalResolution(ScenePointCloud, SceneNeighbors.Connect);

ScenePointCloud = findMeshNormals(ScenePointCloud);

ScenePointCloud.Normal( ScenePointCloud.Normal * [0;0;1] > 0,:) = - ScenePointCloud.Normal(ScenePointCloud.Normal * [0;0;1]> 0,:);


options.hs = ScenePointCloud.Resolution/2;

tau = options.hs^2/4;

ItL = makeExplicitLaplaceBeltrami( strcat(ProjectRoot,'/models/object/itokawa/meshlab/Mesh3D.off') , options, BDF, tau, 1);


SceneScaleParameter = findScaleParamter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

ScaleParameterAbsolute = bsxfun(@plus, SceneScaleParameter, ScenePointCloud.Resolution);

SceneSignal = performBDFDiffusion(ScenePointCloud.Signal, NumSteps, ItL);

SceneDoG = buildDoG(SceneSignal, SceneScaleParameter, DoGNormalize);

SceneKeypoint = findKeypoint(SceneDoG, ScenePointCloud, SceneScaleParameter, SceneNeighbors.Distance, KeypointMethod, CompareMethod);

SceneKeypoint.Scale = SceneKeypoint.Scale + ScenePointCloud.ResolutionLocal(SceneKeypoint.LocationIndex);


SceneNMSKeypoint = applyNMS(ScenePointCloud, SceneDoG, SceneKeypoint, t_scale, t_range, DoGNormalize, CompareMethod);
        
SceneFeatures = buildHistogram(ScenePointCloud, SceneNMSKeypoint);


        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Match features from scene to model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SceneFeatures  ModelFeatures

[a,b]=max(SceneFeatures.Scale);

figure
for i = 1 : SceneFeatures.Count
    clf
    plot([SceneFeatures.HistAlpha(i,:), SceneFeatures.HistPhi(i,:), SceneFeatures.HistBeta(i,:), SceneFeatures.HistTheta(i,:)],'linewidth',4)
    xlabel(sprintf('Features $\\alpha, \\phi, \\beta, \\theta$. 45 bins each'),'interpreter', 'latex')
    ylabel('Percentage')
    title('Example Feature Descriptor')
    ax = gca;
    xlim([1 180])
    ax.XTick = [1, 45, 90, 135, 180];
    ax.YAxis.FontSize = 50;
    ax.XAxis.FontSize = 50;
    drawnow
    pause(0.5)
end



% 
% ScenePointsInModelFrame = bsxfun(@minus, (Rot' * ScenePointCloud.Location')', (Rot' * trans')');
% 
% ClosestModel2ScenePoints = rangesearch(ModelPointCloud.Location, ScenePointsInModelFrame, ScenePointCloud.Resolution);
% 
% c = 0;
% for i = 1 : length(ClosestModel2ScenePoints)
%     if ~isempty(ClosestModel2ScenePoints{i,1})
%         c = c + 1;
%         TruePointCorr(c, 1:2) = [i, ClosestModel2ScenePoints{i,1}(1)];
%     end
% end


t_scale = 0.5;
t_dist = 0.2;
Correspondences = matchDescriptors(SceneFeatures, ModelFeatures, t_scale, t_dist)



figure
plot3(ModelPointCloud.Location(:,1),ModelPointCloud.Location(:,2),ModelPointCloud.Location(:,3),'r.')
hold on
plot3(ScenePointCloud.Location(:,1),ScenePointCloud.Location(:,2),ScenePointCloud.Location(:,3),'b.')
hold on
plot3(SceneFeatures.Location(:,1),SceneFeatures.Location(:,2),SceneFeatures.Location(:,3),'r.','markersize',20)
axis off
axis equal
view(0,-67)
plot3(ModelFeatures.Location(Correspondences.Features(:,2),1),...
      ModelFeatures.Location(Correspondences.Features(:,2),2),...
      ModelFeatures.Location(Correspondences.Features(:,2),3),'b.','markersize',20)
for i = 1 :1000 : Correspondences.Count
    plot3([SceneFeatures.Location(Correspondences.Features(i,1),1), ModelFeatures.Location(Correspondences.Features(i,2),1)],...
          [SceneFeatures.Location(Correspondences.Features(i,1),2), ModelFeatures.Location(Correspondences.Features(i,2),2)],...
          [SceneFeatures.Location(Correspondences.Features(i,1),3), ModelFeatures.Location(Correspondences.Features(i,2),3)],...
          'g-')
    
%     plot3([ScenePointCloud.Location(Correspondences.Features(i,1),1), ModelPointCloud.Location(Correspondences.Features(i,2),1)],...
%           [ScenePointCloud.Location(Correspondences.Features(i,1),2), ModelPointCloud.Location(Correspondences.Features(i,2),2)],...
%           [ScenePointCloud.Location(Correspondences.Features(i,1),3), ModelPointCloud.Location(Correspondences.Features(i,2),3)],...
%           'g-')
end



% 
% [Pose, BestInlierRatio, BestCorrespond, BestPose] = matchDescriptorsRANSAC(Correspondences, SceneFeatures, ModelFeatures, ScenePointCloud, ModelPointCloud);
% ScenePointsInModelFrame = bsxfun(@minus, (BestPose.Rotation' * ScenePointCloud.Location')', (BestPose.Rotation' * BestPose.translation')');
% figure
% plot3(ModelPointCloud.Location(:,1),ModelPointCloud.Location(:,2),ModelPointCloud.Location(:,3),'r.')
% hold on
% plot3(ScenePointsInModelFrame(:,1),ScenePointsInModelFrame(:,2),ScenePointsInModelFrame(:,3),'b.')
% axis off
% axis equal
% view(0,-67)


% figure
% plot3(ScenePointCloud.Location(:,1),ScenePointCloud.Location(:,2),ScenePointCloud.Location(:,3),'r.')
% hold on
% plot3(ScenePointCloud.Location(488,1),ScenePointCloud.Location(488,2),ScenePointCloud.Location(488,3),'b.','markersize',20)




[GeometricCorrespondence, GroupCorrespondence, GroupValues] = groupKeypoints(Correspondences, SceneFeatures, ModelFeatures, t_gc, t_corr);



GroupPose = findRigidPose(GroupCorrespondence, GroupValues, ModelFeatures, SceneFeatures);


clear PoseEstError
for i =  1: length(GroupPose)

    GroupPose{i,1}.Rotation
    GroupPose{i,1}.translation
    [Pose1, movingReg, rmse] = pcregrigid(pointCloud(ModelPointCloud.Location),...
                               pointCloud(ScenePointCloud.Location),'MaxIter', 30, 'InitialTransform',...
                               affine3d([GroupPose{i,1}.Rotation,[0;0;0]; GroupPose{i,1}.translation, 1]));
                           rmse

%     PoseEstError(i,1:3) = GroupPose{i,1}.translation - trans;
%     GroupPose{i,1}.Rotation
    ScenePointsInModelFrameUpdate = bsxfun(@minus, (Pose1.T(1:3,1:3)' * ScenePointCloud.Location')', (Pose1.T(1:3,1:3)' * Pose1.T(4,1:3)')');
i
    ScenePointsInModelFrame = bsxfun(@minus, (GroupPose{i,1}.Rotation' * ScenePointCloud.Location')', (GroupPose{i,1}.Rotation' * GroupPose{i,1}.translation')');
%     ScenePointsInModelFrame = bsxfun(@minus, (GroupPose{i,1}.Rotation * ScenePointCloud.Location')', -GroupPose{i,1}.translation);
    ScenePointsInModelFrameTrue = bsxfun(@minus, (Rot' * ScenePointCloud.Location')', (Rot' * trans')');
%     
    clf
    subplot(1,2,1)
    plot3(ModelPointCloud.Location(:,1),ModelPointCloud.Location(:,2),ModelPointCloud.Location(:,3),'r.')
    hold on
    plot3(ScenePointsInModelFrame(:,1),ScenePointsInModelFrame(:,2),ScenePointsInModelFrame(:,3),'b.')
    axis off
    axis equal
    view(-90,0)
    drawnow
    
    subplot(1,2,2)
    plot3(ModelPointCloud.Location(:,1),ModelPointCloud.Location(:,2),ModelPointCloud.Location(:,3),'r.')
    hold on
%     plot3(ScenePointsInModelFrameTrue(:,1),ScenePointsInModelFrameTrue(:,2),ScenePointsInModelFrameTrue(:,3),'b.')
    plot3(ScenePointsInModelFrameUpdate(:,1),ScenePointsInModelFrameUpdate(:,2),ScenePointsInModelFrameUpdate(:,3),'b.')
    axis off
    axis equal
    view(-90,0)
    drawnow
    pause
    
end


% 
% 
% figure
% plot3(ModelPointCloud.Location(:,1),ModelPointCloud.Location(:,2),ModelPointCloud.Location(:,3),'r.')
% hold on
% quiver3(ModelPointCloud.Location(:,1),ModelPointCloud.Location(:,2),ModelPointCloud.Location(:,3),ModelPointCloud.Normal(:,1),ModelPointCloud.Normal(:,2),ModelPointCloud.Normal(:,3),10)
% axis off
% quiver3(ScenePointsInModelFrame(:,1),ScenePointsInModelFrame(:,2),ScenePointsInModelFrame(:,3),ScenePointCloud.Normal(:,1),ScenePointCloud.Normal(:,2),ScenePointCloud.Normal(:,3),10)


% 
% 
% figure
% plot3(ModelPointCloud.Location(:,1),ModelPointCloud.Location(:,2),ModelPointCloud.Location(:,3),'r.')
% axis off
% hold on
% plot3(ModelNMSKeypoint.Location(:,1),ModelNMSKeypoint.Location(:,2),ModelNMSKeypoint.Location(:,3),'b.','markersize',20)
% ScenePointsInModelFrame = bsxfun(@minus, (Rot' * ScenePointCloud.Location')', (Rot' * trans')');
% plot3(ScenePointsInModelFrame(:,1),ScenePointsInModelFrame(:,2),ScenePointsInModelFrame(:,3),'y.')
% SceneKeypointsInModelFrame = bsxfun(@minus, (Rot' * SceneNMSKeypoint.Location')', (Rot' * trans')');
% 
% plot3(SceneKeypointsInModelFrame(:,1),SceneKeypointsInModelFrame(:,2),SceneKeypointsInModelFrame(:,3),'g.','markersize',20)
% 
% 
% 




























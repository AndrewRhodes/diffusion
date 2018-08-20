
close all
clear
clc

global ProjectRoot; % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


alpha = 1;

porder = 4; % order of interpolation
dim = 3; % dimension
Lorder = 2; % Cartesian Laplace order
spacing = 0.1; % spacing of embedding grid


ShowPlot = 0;
Model = 'armadillo/Armadillo_e1_100000';

BDF = 4;
tauFraction = 1/8;
maxTauNumer = 2.5;
DoGNormalize = 1;

% NumIter = 50;
% NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];

t_scale = 0.7;
t_DoG = 0.9;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(Model,'.ply');
FileNameModelOff = strcat(Model,'.off');


load('Curvature.mat')
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
NormalRotations = findNormalsRotation(PointCloud.Normal);



% % % % % % % % % %
tau = spacing * tauFraction;
MaxTau = maxTauNumer / spacing;
NumSteps = round(MaxTau);
% % % % % % % % % %


[Neighbors, NeighborFaces] = findAdjacentNeighbors(PointCloud);

% [PK1, PK2, PD1, PD2, MK, GK] = findPointCurvatures(PointCloud, NormalRotations, Neighbors.Connect);

PointCloud.Signal = MK;
% clear PK1 PK2 PD1 PD2 GK NeighborFaces
% stdMK = std(MK);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                                      
[ItL, Eplot, CP, CPFACE] = makeImplicitLaplaceBeltrami( PointCloud, porder, Lorder, spacing, dim, BDF, tau, alpha);

FaceInterpolateWeights = interpBarycenterTriangle(PointCloud, CP, CPFACE);

% CPSignal = FaceInterpolateWeights * PointCloud.Signal;

ScaleParameter = findScaleParamter(tau, alpha, NumSteps, 'Laplacian', 'Natural');


save ItL ItL
save Eplot Eplot
sace CP CP
save CPFACE CPFACE
save FaceInterpolateWeights FaceInterpolateWeights


% CPSignal = FaceInterpolateWeights * PointCloud.Signal;
% 
% CPSignalDiffused = performBDFDiffusion(CPSignal, NumSteps, ItL);
% 
% Signal = performEplotProjection(CPSignalDiffused, Eplot);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Difference of Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect Extrema
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% [Keypoint, DoGMaximum] = findKeypointDoGMax(DoG, ScaleParameter, PointCloud, Neighbors.Distance, 'Old');
        
        
        


% vidObj = VideoWriter('DIArmadilloVideoMeshLab.avi');
% vidObj.FrameRate = 90;
% open(vidObj);
% % [SphereX, SphereY, SphereZ] = sphere(100);
% % 
% % SphereX = reshape(SphereX,[],1);
% % SphereY = reshape(SphereY,[],1);
% % SphereZ = reshape(SphereZ,[],1);
% Min = 200;
% for i = 1 : tauNumerator-1
% %     i = i +1;
% %      hh = surfl(reshape(SphereAtPoint(:,PlotOrder(1)),101,101), reshape(SphereAtPoint(:,PlotOrder(2)),101,101), reshape(SphereAtPoint(:,PlotOrder(3)),101,101));
% %     set(hh, 'FaceColor','b', 'EdgeColor','none', 'FaceAlpha',0.5)
%     i
%     figure('Units', 'Normalized', 'OuterPosition', [0 0 0.98 0.98])
% %     colormap('gray')
%     trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), Signal(:,i), 'EdgeColor', 'none');
%     view(180, -90)
% 	ax = gca;
%     axis equal
%     axis off
%     hold on
%     CurrentDoGMax = DoGMaximum.LocationCell{i,1};
%     plot3(PointCloud.Location(CurrentDoGMax,1), PointCloud.Location(CurrentDoGMax,2),...
%         PointCloud.Location(CurrentDoGMax,3), 'k.','MarkerSize', 10)
%     if (i > Min) 
%         CurrentKeypoints = vertcat(Keypoint.LocationCell{Min:i,1});
% %         for j = 1 : length(CurrentKeypoints)
% %             SphereAtPoint = bsxfun(@plus, Keypoint.Scale(Keypoint.Location*[SphereX, SphereY, SphereZ], PointCloud.Location(FeaturePoint.Location(SortOrder(i)), :));
%         plot3(PointCloud.Location(CurrentKeypoints,1), PointCloud.Location(CurrentKeypoints,2),...
%             PointCloud.Location(CurrentKeypoints,3), 'r.','Markersize',50)
% %         end
%     elseif i < Min
%         CurrentKeypoints = vertcat(Keypoint.LocationCell{i,1});
% %         for j = 1 : length(CurrentKeypoints)
% %             SphereAtPoint = bsxfun(@plus, Keypoint.Scale(Keypoint.Location*[SphereX, SphereY, SphereZ], PointCloud.Location(FeaturePoint.Location(SortOrder(i)), :));
%         plot3(PointCloud.Location(CurrentKeypoints,1), PointCloud.Location(CurrentKeypoints,2),...
%             PointCloud.Location(CurrentKeypoints,3), 'r.','Markersize',50)
%     end
%     
%     cbar = colorbar('FontSize',15);
%     caxis([-0.15,0.2])
%     ylabel(cbar, 'Mean Curvature', 'FontSize',30);
%     ax.XLim = [-240 240];
%     ax.YLim = [-115 200];
%     drawnow
%     pause(0.5)
%     
%     FileLocation = strcat(ProjectRoot,'/main/DI/keypointdata/armadillo/Images/First/ArmadilloStep_',num2str(i));
%     savefig(gcf, FileLocation)
%         
%     writeVideo(vidObj, getframe(gca));
%     
%     close(gcf);
%     
%     
%     
% end
% close(vidObj)











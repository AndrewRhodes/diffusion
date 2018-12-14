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
options.dtype = 'geodesic';
% options.dtype = 'euclidean';

ShowPlot = 0;
ModelFolder = 'armadillo/';
Model = 'Armadillo_e1_100000';

BDF = 2;
tauFraction = 1/10;
tauNumerator = 3000;

KeypointSearchMethod = 'Global'; % 'Local', 'Global'
Radius = 1;


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

[DoG, AbsDoG] = buildDoGNew(Signal, ScaleParameter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect Extrema
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [Keypoint, Zeropoint] = findKeypointNew(PointCloud, ScaleParameter, Neighbors.Connect, DoG, AbsDoG);
%         
% NewKeypoint = checkKeypointSphere(PointCloud, ScaleParameter, Keypoint, Zeropoint, Radius, KeypointSearchMethod);


% save(strcat(Model,'_N',num2str(NumSteps),'_Keypoint.mat'),'Keypoint','-v7.3')
% save(strcat(Model,'_N',num2str(NumSteps),'_NewKeypoint.mat'),'NewKeypoint','-v7.3')
% save(strcat(Model,'_N',num2str(NumSteps),'_Zeropoint.mat'),'Zeropoint','-v7.3')


load(strcat(Model,'_N',num2str(NumSteps),'_Keypoint.mat'),'Keypoint')
load(strcat(Model,'_N',num2str(NumSteps),'_NewKeypoint.mat'),'NewKeypoint')
load(strcat(Model,'_N',num2str(NumSteps),'_Zeropoint.mat'),'Zeropoint')


% ZeroLogic = (Keypoint.Scale ./ PointCloud.ResolutionLocal(Keypoint.Location)) < 4;
ZeroLogic = (Keypoint.Scale - PointCloud.ResolutionLocal(Keypoint.Location)) ./ PointCloud.Resolution < 3;
Keypoint.Scale(ZeroLogic) = [];
Keypoint.Location(ZeroLogic) = [];
Keypoint.Level(ZeroLogic) = [];
Keypoint.Sign(ZeroLogic) = [];
Keypoint.Count = length(Keypoint.Scale);


% ZeroLogic = (NewKeypoint.Scale ./ PointCloud.ResolutionLocal(NewKeypoint.Location))  < 4;
ZeroLogic = (NewKeypoint.Scale - PointCloud.ResolutionLocal(NewKeypoint.Location)) ./ PointCloud.Resolution < 3;
NewKeypoint.Scale(ZeroLogic) = [];
NewKeypoint.Location(ZeroLogic) = [];
NewKeypoint.Level(ZeroLogic) = [];
NewKeypoint.Sign(ZeroLogic) = [];
NewKeypoint.Count = length(NewKeypoint.Scale);

UniqueLevelNewKeypoint = unique(NewKeypoint.Level);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the image for the video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sunvector = [0;0;-1];
% sunvector = sunvector / norm(sunvector);
% intensity = 0.9*PointCloud.Normal*sunvector;
% intensity(intensity<0) = 0;


clear vidObj
vidObj = VideoWriter('Armadillo_NewKeypoints.avi');
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

    i=586
%     subplot(1,2,1)
    t1 = trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), Signal(:,i), 'EdgeColor', 'none');
%     t1 = trisurf(PointCloud.Face, PointCloud.Location(:,3), PointCloud.Location(:,1), PointCloud.Location(:,2), Signal(:,i), 'EdgeColor', 'none');    
    ax = gca;
    view(180,-90)
%     view(ax,-130,0) 
%     view(ax,-90,0) % 
    axis equal
    axis off
    hold on
    FigSize1 = ax.Position;
    cbar = colorbar('FontSize',15);
    
    
%     subplot(1,2,2)
%     t2 = trisurf(PointCloud.Face, PointCloud.Location(:,3), PointCloud.Location(:,1), PointCloud.Location(:,2), Signal(:,i), 'EdgeColor', 'none');    
%     ax2 = gca;
%     view(ax2,-60, 0)
%     axis equal
%     axis off
%     hold on
%     FigSize2 = ax2.Position;
%     cbar = colorbar('FontSize',15);
%     ax2.Position = [FigSize2(1:2), FigSize1(3:4)];

    
    for i = 586 %i = [105,121,135,141,142,144,161,164,169,171,236,240,344,365,367,533,586,636,747,843,1516,1652,1895,2254,2605] % 
         
    CurrentNewKeypoints = NewKeypoint.Location( NewKeypoint.Level == i);
    CurrentNewScales = NewKeypoint.Scale(NewKeypoint.Level == i);
    CurrentSigns = NewKeypoint.Sign(NewKeypoint.Level == i);
    for j = 1 : length(CurrentNewKeypoints)
        j
        
    
        SphereAtPoint = bsxfun(@plus, CurrentNewScales(j)*[SphereX, SphereY, SphereZ], PointCloud.Location(CurrentNewKeypoints(j),:));
        hh = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
        set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)
        
        if CurrentSigns(j) < 0
            plot3(PointCloud.Location(CurrentNewKeypoints,1), PointCloud.Location(CurrentNewKeypoints,2),...
        PointCloud.Location(CurrentNewKeypoints,3), 'b.','MarkerSize', 30);
%             set(hh, 'FaceColor','b', 'EdgeColor','none', 'FaceAlpha',0.5)
        else
            plot3(PointCloud.Location(CurrentNewKeypoints,3), PointCloud.Location(CurrentNewKeypoints,1),...
        PointCloud.Location(CurrentNewKeypoints,2), 'r.','MarkerSize', 30);
%             set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)
        end
        drawnow
    end
    end

    
%     CurrentZeropoints = Zeropoint.Location( Zeropoint.Level == i);
%     p1 = plot3(PointCloud.Location(CurrentZeropoints,1), PointCloud.Location(CurrentZeropoints,2),...
%         PointCloud.Location(CurrentZeropoints,3), 'k.','MarkerSize', 10);
    
    
    
    CurrentKeypointsLow = Keypoint.Location( ((Keypoint.Level == i) & (Keypoint.Sign < 0)) );
    p2 = plot3(PointCloud.Location(CurrentKeypointsLow,1), PointCloud.Location(CurrentKeypointsLow,2),...
        PointCloud.Location(CurrentKeypointsLow,3), 'b.','MarkerSize', 10);
    
    
    
    CurrentKeypointsHigh = Keypoint.Location( ((Keypoint.Level == i) & (Keypoint.Sign > 0)) );
    p3 = plot3(PointCloud.Location(CurrentKeypointsHigh,1), PointCloud.Location(CurrentKeypointsHigh,2),...
        PointCloud.Location(CurrentKeypointsHigh,3), 'r.','MarkerSize', 10);
    
    ylabel(cbar, 'Mean Curvature', 'FontSize',30);
    ax.XLim = [-240 250];
    ax.YLim = [-150 250];
    ax.ZLim = [-200, 300];
    title(sprintf('Diffusion Step: %d, Absolute Scale: %0.2f',i, ScaleParameterAbsolute(i+1)),'fontsize',25)
    caxis([-0.15,0.2])
%     if ~any((i == 1) || (i == tauNumerator-1))
%         if isempty(CurrentNewKeypoints)
%             [~, hobj, ~, ~] = legend([t1, p1, p2, p3],{'Model','Edge','Min','Max'},'Fontsize',35,'Location','West');
%             hobj(7).MarkerSize = 80;
%             hobj(9).MarkerSize = 80;
%             hobj(11).MarkerSize = 80;
%         else
%             [~, hobj, ~, ~] = legend([t1, p1, p2, p3, hh],{'Model','Edge','Min','Max', 'Keypoint'},'Fontsize',35,'Location','West');
%             hobj(8).MarkerSize = 80;
%             hobj(10).MarkerSize = 80;
%             hobj(12).MarkerSize = 80;
%         end
%     end
    
    drawnow
%     pause(0.1)
    
%     FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'Video/Images/NewRun/Step_',num2str(i));
%     savefig(fig, FileLocation)
        
%     writeVideo(vidObj, getframe(fig));
    
    clf(fig);
    
end
close(fig);

close(vidObj)




    






















    

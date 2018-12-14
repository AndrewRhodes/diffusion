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
ModelFolder = 'bunny/';
Model = 'Bunny_e1';

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
vidObj = VideoWriter('Bunny_NewKeypoints.avi');
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

    t1 = trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), Signal(:,i), 'EdgeColor', 'none');
%     t1 = trisurf(PointCloud.Face, PointCloud.Location(:,3), PointCloud.Location(:,1), PointCloud.Location(:,2), Signal(:,i), 'EdgeColor', 'none');    
    ax = gca;
    view(0, 90)
    axis equal
    axis off
    hold on
    FigSize1 = ax.Position;
    cbar = colorbar('FontSize',15);
    
         
    CurrentNewKeypoints = NewKeypoint.Location( NewKeypoint.Level == i);
    CurrentNewScales = NewKeypoint.Scale(NewKeypoint.Level == i);
    CurrentSigns = NewKeypoint.Sign(NewKeypoint.Level == i);
    
    for j = 1 : length(CurrentNewKeypoints)
        SphereAtPoint = bsxfun(@plus, CurrentNewScales(j)*[SphereX, SphereY, SphereZ], PointCloud.Location(CurrentNewKeypoints(j),:));
        hh = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
        set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)
        
        if CurrentSigns(j) < 0
            set(hh, 'FaceColor','b', 'EdgeColor','none', 'FaceAlpha',0.5)
        else
            set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)
        end
    end

    
    CurrentZeropoints = Zeropoint.Location( Zeropoint.Level == i);
    p1 = plot3(PointCloud.Location(CurrentZeropoints,1), PointCloud.Location(CurrentZeropoints,2),...
        PointCloud.Location(CurrentZeropoints,3), 'k.','MarkerSize', 10);
    
    
    
    CurrentKeypointsLow = Keypoint.Location( ((Keypoint.Level == i) & (Keypoint.Sign < 0)) );
    p2 = plot3(PointCloud.Location(CurrentKeypointsLow,1), PointCloud.Location(CurrentKeypointsLow,2),...
        PointCloud.Location(CurrentKeypointsLow,3), 'b.','MarkerSize', 10);
    
    
    
    CurrentKeypointsHigh = Keypoint.Location( ((Keypoint.Level == i) & (Keypoint.Sign > 0)) );
    p3 = plot3(PointCloud.Location(CurrentKeypointsHigh,1), PointCloud.Location(CurrentKeypointsHigh,2),...
        PointCloud.Location(CurrentKeypointsHigh,3), 'r.','MarkerSize', 10);
    
    ylabel(cbar, 'Mean Curvature', 'FontSize',30);
    ax.XLim = [-90 70];
    ax.YLim = [15 135];
    ax.ZLim = [-50 50];
    title(sprintf('Diffusion Step: %d, Absolute Scale: %0.2f',i, ScaleParameterAbsolute(i+1)),'fontsize',25)
    caxis([-0.1202, 0.1793])
    
    drawnow
    
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'Video/Images/NewRun/Step_',num2str(i));
    savefig(fig, FileLocation)
        
    writeVideo(vidObj, getframe(fig));
    
    clf(fig);
    
end
close(fig);

close(vidObj)




    






















    

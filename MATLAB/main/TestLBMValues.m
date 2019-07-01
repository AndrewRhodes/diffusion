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
options.dtype = 'geodesic'; % 'euclidean', 'geodesic' %
options.htype = 'ddr'; % 'psp', 'ddr'

% Destination = 'Mesh_rho4_ddr_geo_sph' %'Mesh_rho4_ddr_geo_t0.25'
ModelFolder = 'bunny/';
Model = 'Bunny_e1';

BDF = 1;
NumIter = 5;
NumSteps = 2000;
DoGNormalize = 'NLoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

t_scale = 1/sqrt(2);
t_range1 = 1/2;
t_range2 = 2;

NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];

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

FileNameItL = strcat(Model,'_ItL','_BDF',num2str(BDF),'_rho',...
              num2str(options.rho),'_',options.dtype(1:3),'_',...
              options.htype,'_t0.25e','_a',num2str(alpha),'.mat'); 
          
setTau = @(e_bar) 0.25*e_bar;
setHs = @(e_bar) 2*e_bar^(1/5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PointCloud.Location, PointCloud.Face, PointCloud.Normal, PointCloud.Signal]...
                = read_ply_all_elements( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
PointCloud = findMeshResolution(PointCloud, 'Model');



FileNameNeighbors = strcat(Model,'_Neighbors.mat');
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
        


% [PK1, PK2, PD1, PD2, MK, GK] = findPointCurvatures(PointCloud, NormalRotations, Neighbors.Connect);
% clear PK1 PK2 PD1 PD2 GK NeighborFaces
Quants = quantile(PointCloud.Signal, [0.25,0.5,0.75]);
MKQuant = PointCloud.Signal;
OutOfBounds = (MKQuant > (Quants(3) + 1.5*(Quants(3)-Quants(1)))) | (MKQuant < (Quants(1) - 1.5*(Quants(3)-Quants(1))));
MKQuant(OutOfBounds) = [];
stdMK = std(MKQuant);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


save_off(PointCloud.Location*4, PointCloud.Face, strcat(FileLocationModel, ModelFolder, Model,'_large','.off') )


tau = setTau(PointCloud.Resolution);

[ItL, LBM, A, Weights] = makeMeshLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), options, BDF, tau, alpha);
ItL = ItL{1};
[ItL_cot, LBM_cot, A_cot, W_cot] = makeCotangentLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), BDF, tau, alpha);
ItL_cot = ItL_cot{1};
[ItL_umb, LBM_umb] = makeUmbrellaLaplaceBeltrami(PointCloud, Neighbors.Connect, tau, alpha, BDF);
ItL_umb = ItL_umb{1};




tau = setTau(PointCloud.Resolution*4);

[ItL_large, LBM_large, A_large, Weights_large] = makeMeshLaplaceBeltrami( strcat(FileLocationModel, ModelFolder, Model,'_large','.off') , options, BDF, tau, alpha);
ItL_large = ItL_large{1};

[ItL_cot_large, LBM_cot_large, A_cot_large, W_cot_large] = makeCotangentLaplaceBeltrami( strcat(FileLocationModel, ModelFolder, Model,'_large','.off'), BDF, tau, alpha);
ItL_cot_large = ItL_cot_large{1};



% if ~exist(strcat(FileLocationMeshItL, FileNameItL), 'file')
%     ItL = makeMeshLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), options, BDF, tau, alpha);
%     save(strcat(FileLocationMeshItL, FileNameItL),'ItL','-v7.3');
% else
%     load(strcat(FileLocationMeshItL, FileNameItL), 'ItL');
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ScaleParameter = findScaleParameter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

GaussianScaleParameter = findScaleParameter(sqrt(2*1*0.25*PointCloud.Resolution), alpha, NumSteps, 'Gaussian', 'Natural');


GaussianScaleParameter ./ (sqrt(2*1*0.25*PointCloud.Resolution*(1:NumSteps)'))


sigma = sqrt(2*1*0.25*PointCloud.Resolution)

t = sigma^2 ./ (2*alpha*(1:NumSteps)')

sqrt(2*1*t)

sqrt(t)

diff(sqrt(2*1*t))

sig = sqrt(2*alpha*tau*(1:NumSteps)')

sig ./ (sqrt(2*1*0.25*PointCloud.Resolution*(1:NumSteps)'))






for i = 1
    
    sprintf('No Noise : %d',i)
    Signal = zeros(PointCloud.LocationCount, NumSteps);
    Signal(:,1) = PointCloud.Signal;
    
    for kk = 1 : NumSteps - 1

        [Signal(:,kk+1), flag] = bicgstab(ItL, Signal(:,kk), 1e-10, 100);

        if flag
            disp(flag)
        end

    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find Difference of Gaussian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detect Extrema
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Distance);
    Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
    Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
    Keypoint.Scale = ScaleParameter(Keypoint.Level);
        
    
    NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range1, 'sigma', DoGNormalize, CompareMethod);
 
 
    
end






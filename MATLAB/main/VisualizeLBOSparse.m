% Andrew Rhodes
% ASEL
% December 2018


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

ModelFolder = 'itokawa/';
Model = 'Itokawa_e1_80000';

% modelFolder = {'bunny','armadillo','buddha','dragon','itokawa'}; % 
% modelNames = {'Bunny_e1','Armadillo_e1_100000','Buddha_e1_50000', 'Dragon_e1_50000', 'Itokawa_e1_80000'};
 


BDF = 1;
NumSteps = 2000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model File Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');
FileNameModelOff = strcat(ModelFolder,Model,'.off');
FileLocationWD = '/media/andrew/WDRhodes/diffusiondata/';

TmpLocation = strcat(ProjectRoot,'/models/object/',ModelFolder,'meshlab/');

FileLocationNeighbors = strcat(FileLocationWD,ModelFolder,'neighbors/');

          
setTau = @(e_bar) 0.25*e_bar;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       

[ItL_mesh, LBM_mesh] = makeMeshLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), options, BDF, tau, alpha);
p_mesh = symrcm(LBM_mesh);


[ItL_Umb, LBM_umb] = makeUmbrellaLaplaceBeltrami(PointCloud, Neighbors.Connect, tau, alpha, BDF);
p_umb = symrcm(LBM_umb);


% [ItL_cot, LBM_cot] = makeCotangentLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), BDF, tau, alpha);
% p_cot = symrcm(LBM_cot);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures

figure
spy(LBM_mesh(p_mesh,p_mesh))
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'xlabel',[])

% figure
% spy(LBM_umb(p_umb,p_umb))
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
% set(gca,'xlabel',[])


% figure
% spy(LBM_cot(p_cot,p_cot))
% axis off







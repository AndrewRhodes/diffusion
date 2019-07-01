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

Destination = 'Mesh_rho4_ddr_geo' %'Mesh_rho4_ddr_geo_t0.25'
ModelFolder = 'armadillo/';
Model = 'Armadillo_e1_100000';

NumIter = 15;
DoGNormalize = 'DoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'

l_range1 = 1/2;
l_range2 = 2;
k = 1.2;%2^(1/4);
l_scale = 1/k;
l_ebar = 80;


set_ltau = @(e_bar) 4*e_bar;
set_hs = @(e_bar) 2*e_bar^(1/5);
set_t0 = @(e_bar) e_bar/4 ;
set_NumSteps = @(e_bar, t_0) ceil(log((l_ebar*e_bar)^2 / (2*alpha*t_0)) / (2*log(k)));


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

FileNameLBM = strcat(Model,'_LBM','_mesh_', options.dtype(1:3),'_rho',...
              num2str(options.rho),'_', options.htype,'.mat'); 
          


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [PointCloud.Location, PointCloud.Face, PointCloud.Normal, PointCloud.Signal]...
%                 = read_ply_all_elements( fullfile( FileLocationModel, strcat(ModelFolder,'Armadillo_e1_100000_GK2','.ply') ) );
            
% %  PointCloud.Color = .2126 * (PointCloud.Color(:,1)/255).^2.2 + .7152 * (PointCloud.Color(:,2)/255).^2.2 + .0722 * (PointCloud.Color(:,3)/255).^2.2;
% %            
% % PointCloud.Signal = PointCloud.Color;

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
% tau = set_tau(PointCloud.Resolution);
l_tau = set_ltau(PointCloud.Resolution);
t_0 = set_t0(PointCloud.Resolution);
NumSteps = set_NumSteps(PointCloud.Resolution, t_0);

if strcmp(options.htype, 'psp')
    options.hs = set_hs(PointCloud.Resolution);
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


if ~exist(strcat(FileLocationMeshItL, FileNameLBM), 'file')
    [LBM] = makeMeshLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), options);
    save(strcat(FileLocationMeshItL, FileNameLBM), 'LBM','-v7.3');
else
    load(strcat(FileLocationMeshItL, FileNameLBM), 'LBM');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[tn, tau_n, sigma_n] = findScaleStep(k, t_0, alpha, NumSteps);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusion of Mean Curvature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : length(NoiseVec)
	
    for i = 1 : NumIter
        
        sprintf('Std %0.1f : %d',NoiseVec(j),i)
        
        SignalNoisy = PointCloud.Signal + NoiseVec(j)*stdMK*randn(PointCloud.LocationCount,1);
   
        [Signal, IterCount] = performBDFDiffusion(SignalNoisy, LBM, alpha, tau_n, NumSteps, l_tau);
        
        
%         Signal = performBDFDiffusion(SignalNoisy, NumSteps, ItL);

%         Signal = performBDFDiffusion_cpp(ItL, SignalNoisy, NumSteps);
  
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find Difference of Gaussian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        DoG = buildDoG(Signal, sigma_n, DoGNormalize);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Detect Extrema
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Distance);
        Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
        Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
        Keypoint.Scale = sigma_n(Keypoint.Level);
        Keypoint.Count = length(Keypoint.Scale);
        
        

%         ZeroLogic = Keypoint.Scale < 1* PointCloud.Resolution;
% 
%         FNames = fieldnames(Keypoint);
%         for jj = 1 : length(FNames)
%             if strcmpi(FNames{jj},'Count')
%                 Keypoint = rmfield(Keypoint, FNames{jj});
%             else
%                 Keypoint.(FNames{jj})(ZeroLogic,:) = [];
%             end
%         end
%         Keypoint.Count = length(Keypoint.LocationIndex)



%         LevelTable = tabulate(Keypoint.Level);
%         figure
%         plot(sigma_n(LevelTable(:,1)), ...
%              LevelTable(:,2), 'k-', 'Linewidth', 4)

 
 
        NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, l_scale, l_range1, 'sigma', DoGNormalize, CompareMethod);               
        
        % SubKeypoint = findSubKeypoint(Keypoint, ScaleParameterAbsolute, DoG, PointCloud, Neighbors, NeighborFaces);
        
        
        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'SignalNoise/',Destination,'/Std_',num2str(NoiseVec(j)),'/');
        FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
        
        
        FileName = strcat('NMSKeypoint','_sigma_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
        
        
        NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, l_scale, l_range2, 'ebar', DoGNormalize, CompareMethod);
        FileName = strcat('NMSKeypoint','_ebar_Iter',num2str(i),'.mat');
        save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
    end
end




for i = 1
    
    sprintf('No Noise : %d',i)
    
%     MK = mean(LBMsym * PointCloud.Location,2);


%     [theta, phi, r] = cart2sph(PointCloud.Location(:,1),...
%         PointCloud.Location(:,2),PointCloud.Location(:,3));
    
    [Signal, IterCount] = performBDFDiffusion(PointCloud.Signal, LBM, alpha, tau_n, NumSteps, l_tau);
    
    
    
%     Signal = performBDFDiffusion(phi, NumSteps, ItL);
%     Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
%     Signal = performBDFDiffusion_cpp(ItL, PointCloud.Signal, NumSteps);
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find Difference of Gaussian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DoG = buildDoG(Signal, sigma_n, DoGNormalize);
%     for d = 1 : PointCloud.LocationCount
%        DoG(d,:) = DoG(d,:) / max(abs(DoG(d,:)));
%     end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detect Extrema
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Distance);
    Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
    Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
    Keypoint.Scale = sigma_n(Keypoint.Level);
    Keypoint.Count = length(Keypoint.Scale);
    
    NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, l_scale, l_range1, 'sigma', DoGNormalize, CompareMethod);
 
    %     SubKeypoint = findSubKeypoint(Keypoint, ScaleParameterAbsolute, DoG, PointCloud, Neighbors.Connect, NeighborFaces.Connect);
    
    
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',ModelFolder,'SignalNoise/',Destination,'/');
    FileName = strcat('Keypoint','.mat');
    
    save(fullfile(FileLocation, FileName), 'Keypoint', '-v7.3')
    
    FileName = strcat('NMSKeypoint_sigma','.mat');
    save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')

    
    NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, l_scale, l_range2, 'ebar', DoGNormalize, CompareMethod);
    FileName = strcat('NMSKeypoint_ebar','.mat');
    save(fullfile(FileLocation, FileName), 'NMSKeypoint', '-v7.3')
  
end


% LevelTable = tabulate(Keypoint.Level);
% figure
% plot(sigma_n(LevelTable(:,1))/PointCloud.Resolution, ...
%      LevelTable(:,2), 'k-', 'Linewidth', 4)
% 
%  
% ZeroLogic = Keypoint.Scale < 2 * PointCloud.Resolution;
% % ZeroLogic = Keypoint.Scale > sigma_n(end-3);
% 
% FNames = fieldnames(Keypoint);
% for jj = 1 : length(FNames)
%     if strcmpi(FNames{jj},'Count')
%         Keypoint = rmfield(Keypoint, FNames{jj});
%     else
%         Keypoint.(FNames{jj})(ZeroLogic,:) = [];
%     end
% end
% Keypoint.Count = length(Keypoint.LocationIndex)
% 
% 
% % figure
% % hold on
% % for i = 1 : Keypoint.Count
% %     plot(sigma_n(1:end-2), DoG(Keypoint.LocationIndex(i),1:end-1),'-')
% % %     pause
% % end
% 
% sunvector = [0;0;-1];
% sunvector = sunvector / norm(sunvector);
% intensity = 0.9*PointCloud.Normal*sunvector;
% intensity(intensity<0) = 0;
% 
% % 
% figure
% colormap('gray')
% % trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), 'EdgeColor', 'none');
% trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), intensity, 'EdgeColor', 'none');
% view(180, -90)
% axis equal
% axis off
% hold on
% 
% 
% 
% [SphereX, SphereY, SphereZ] = sphere(100);
% SphereX = reshape(SphereX,[],1);
% SphereY = reshape(SphereY,[],1);
% SphereZ = reshape(SphereZ,[],1);
% 
% [Value, LocationIndex] = sort(Keypoint.Scale,'ascend');
% 
% for i = 1:48
% 	SphereAtPoint = bsxfun(@plus, Keypoint.Scale(LocationIndex(i))*[SphereX, SphereY, SphereZ], Keypoint.Location(LocationIndex(i),:));
%     hh = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
%     set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
DE_Armadillo_e1_euc
DE_Armadillo_e1_Umbrella
DE_Armadillo_e1_Cotangent

% 
% 
DE_Armadillo_e1_VNoise
DE_Armadillo_e1_VNoise_euc
DE_Armadillo_e1_VNoise_Umbrella
DE_Armadillo_e1_VNoise_Cotangent

% 
% clear
% 
% DE_Buddha_e1
% DE_Bunny_e1
% DE_Dragon_e1
% DE_Itokawa_e1
% 
% clear



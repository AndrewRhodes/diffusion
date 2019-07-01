% Andrew Rhodes
% ASEL
% Jan. 2019

% Perform scale space on an 2D image of an camera man and find key points.
% Sample the image onto a 3D explicit surface.
% Perform implicit surface diffusion, find key points, and compare.


close all
clear 
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional Paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ProjectRoot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Overall settings
NumSteps = 500;
FileLocation = strcat(ProjectRoot,'/models/image/');
ImageChoice = 'sunflowers' % 'comet67p' %
% FileName = 'Comet67PPlaneImage.off';



% Settings: 2D image diffusion
gausswindow = 5;
sigma2D = 0.75;
DoGNormalize2D = 'DoG';



% Settings: 3D surface diffusion
DoGNormalize3D = 'NLoG';
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'
alpha = 1;
BDF = 1;

t_scale = 1/sqrt(2);
t_range1 = 1/2;
t_range2 = 2;



% Settings: 3D mesh LBO
options.rho = 4;
options.dtype = 'euclidean'; % 'euclidean', 'geodesic' %
options.htype = 'ddr'; % 'psp', 'ddr'





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image Selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(ImageChoice, 'sunflowers')
    
    ImageName = 'sunflowers.png';
    Image = double(rgb2gray(imread(ImageName)));
    Image = Image(:,200:500);
    Image = imresize(Image,0.75);
    Image = Image./255;
   
elseif strcmpi(ImageChoice, 'comet67p')
    
    ImageName = 'Comet67P813x1128.png';
    Image = double(imread(strcat(FileLocation,ImageName)));
    Image = imresize(Image,0.25);
    Image = Image./255;
    
end


% figure
% imshow(Image)


ImageSize = size(Image);
NumVertex = numel(Image);
NumPixels = numel(Image);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



ImageFileName = strcat(ImageChoice ,'_e1_', num2str(NumPixels));
ImageFileName3 = strcat(ImageChoice ,'_e3_', num2str(NumPixels));
ImageFileNameOff = strcat(FileLocation, ImageFileName, '.off');
ImageFileNameOff3 = strcat(FileLocation, ImageFileName3, '.off');












%% Scale space on 3D surface


NumPointSurf = length(1:NumVertex);
MiddleSurfPoint = round(NumPointSurf/2);

[xSurf3D, ySurf3D, zSurf3D] = meshgrid(1:ImageSize(2),1:ImageSize(1),0);

PointCloud.Location = [xSurf3D(:), ySurf3D(:), zSurf3D(:)];
PointCloud.LocationCount = length(PointCloud.Location);
PointCloud.Face = delaunay(ySurf3D(:), xSurf3D(:));
PointCloud.FaceCount = length(PointCloud.Face);
PointCloud.FaceArea = findFaceArea(PointCloud.Location, PointCloud.Face);
PointCloud.Signal = reshape(Image, [], 1);
PointCloud = findMeshResolution(PointCloud, 'Model');


FileNameNeighbors = strcat(ImageFileName,'_Neighbors.mat');
if ~exist( strcat( FileLocation, FileNameNeighbors), 'file')   
    [Neighbors, ~, PointCloud] = findAdjacentNeighbors(PointCloud);    
    save(strcat( FileLocation, FileNameNeighbors) ,'Neighbors', '-v7.3')
    Neighbors3D = Neighbors;
    clear Neighbors
else
    load(strcat( FileLocation, FileNameNeighbors), 'Neighbors');
    Neighbors3D = Neighbors;
    clear Neighbors
end

PointCloud = findLocalResolution(PointCloud, Neighbors3D.Distance);


if ~exist(ImageFileNameOff, 'file')
    save_off(PointCloud.Location, PointCloud.Face, ImageFileNameOff);
end

if ~exist(ImageFileNameOff3, 'file')
    save_off(3*PointCloud.Location, PointCloud.Face, ImageFileNameOff3);
end
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Diffusion w/ mesh LBO Euc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setTau = @(e) e;

tau = setTau(PointCloud.Resolution);

ScaleParameter3D = findScaleParameter(tau, 1, NumSteps, 'Laplacian', 'Natural');



[ItL, LBM, A, W] = makeMeshLaplaceBeltrami( ImageFileNameOff, options, BDF, tau, alpha);
ItL_mesh = ItL;
LBM_mesh = LBM;
clear ItL LBM



% Signal_mesh = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL_mesh);
 
Signal_mesh = zeros(PointCloud.LocationCount, NumSteps);
Signal_mesh(:,1) = PointCloud.Signal;
for i = 1 : NumSteps - 1        
    ItL_i = speye(PointCloud.LocationCount, PointCloud.LocationCount) - (1/sqrt(2)^(1/i))*tau * LBM_mesh;
    [Signal_mesh(:,i+1), flag] = bicgstab(ItL_i, Signal_mesh(:,i), 1e-10, 100);    
    if flag
        disp(flag)
    end   
end


% for i = 1 : NumSteps - 1
%     if mod(i,2) % odd
%         if i == 1
%             Sig(:,i) = Signal_mesh(:,i);
%         else
%             Sig(:,i) = Signal_mesh(:,i-1);
%         end
%     else
%         Sig(:,i) = (Signal_mesh(:,i-1) + Signal_mesh(:,i+1))/2;
%     end
% end



DoG_mesh = buildDoG(Signal_mesh, ScaleParameter3D, 'DoG'); % DoGNormalize3D

% for i = 1 : NumSteps - 1
%     imshow(reshape(DoG_mesh(:,i), ImageSize(1), ImageSize(2)),[])
%     drawnow
% end

[Keypoint_mesh.LocationIndex, Keypoint_mesh.Level] = findKeypoint_c(DoG_mesh, Neighbors3D.Distance);
Keypoint_mesh.Location = PointCloud.Location(Keypoint_mesh.LocationIndex,:);
Keypoint_mesh.Scale = ScaleParameter3D(Keypoint_mesh.Level);
Keypoint_mesh.Count = length(Keypoint_mesh.LocationIndex)


LevelTable = tabulate(Keypoint_mesh.Level);
figure
plot(LevelTable(:,1)./PointCloud.Resolution, ...
     LevelTable(:,2)./max(LevelTable(:,2)), 'k-', 'Linewidth', 4)


ZeroLogic = Keypoint_mesh.Scale < 3*PointCloud.Resolution;
         
FNames = fieldnames(Keypoint_mesh);
for jj = 1 : length(FNames)
    if strcmpi(FNames{jj},'Count')
        Keypoint_mesh = rmfield(Keypoint_mesh, FNames{jj});
    else
        Keypoint_mesh.(FNames{jj})(ZeroLogic,:) = [];
    end
end
Keypoint_mesh.Count = length(Keypoint_mesh.LocationIndex)



% NMSKeypoint_mesh = applyNMS(PointCloud, DoG_mesh, Keypoint_mesh, t_scale, t_range1, 'sigma', DoGNormalize3D, CompareMethod);
% figure
% for i = 1 : NumSteps    
%     a = reshape(Signal_mesh(:,i), ImageSize(1), ImageSize(2));
%     imshow(a)    
% end


[SortValue3D, SortOrder3D] = sort(Keypoint_mesh.Scale,'descend');


CircleAngle = linspace(0, 2*pi, 360);
CircleX = cos(CircleAngle);
CircleY = sin(CircleAngle);


PlotFeatures = 1:Keypoint_mesh.Count;%Keypoint_mesh.Count;


figure
imshow(reshape(PointCloud.Signal, ImageSize(1), ImageSize(2)))
hold on
title(sprintf('Mesh: t=%0.2f',tau))

% Still need to work on. 1/3/2019
for i = 1 : length(PlotFeatures)
% % % % %     CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY], [PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),1);PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),2)]);
    
    CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY],...
        [Keypoint_mesh.Location(SortOrder3D(PlotFeatures(i)),2);...
        Keypoint_mesh.Location(SortOrder3D(PlotFeatures(i)),1)]);


    plot(CircleAtPoint(2,:), CircleAtPoint(1,:),'r')
    

%     pause
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% setTau = @(e) 0.5*e;

beta = 3;


PointCloud3 = PointCloud;
PointCloud3.Location = PointCloud3.Location * beta;
PointCloud3 = findLocalResolution(PointCloud3, Neighbors3D.Distance);
PointCloud3 = findMeshResolution(PointCloud3, 'Model')


tau3 = setTau(PointCloud3.Resolution);


% ScaleParameter3D3 = findScaleParameter(beta^2*tau, alpha, NumSteps, 'Laplacian', 'Natural');
ScaleParameter3D3 = findScaleParameter(tau3, alpha, beta*NumSteps, 'Laplacian', 'Natural');

[ItL, LBM, A, W] = makeMeshLaplaceBeltrami( ImageFileNameOff3, options, BDF, tau3, alpha);
ItL_mesh3 = ItL;
LBM_mesh3 = LBM;
clear ItL LBM


% ItL_mesh3{1,1} = speye(PointCloud.LocationCount,PointCloud.LocationCount) - tau3*LBM_mesh3;


Signal_mesh3 = performBDFDiffusion(PointCloud3.Signal, beta*NumSteps, ItL_mesh3);


DoG_mesh3 = buildDoG(Signal_mesh3, ScaleParameter3D3, DoGNormalize3D); % DoGNormalize3D


[Keypoint_mesh3.LocationIndex, Keypoint_mesh3.Level] = findKeypoint_c(DoG_mesh3, Neighbors3D.Distance);
Keypoint_mesh3.Location = PointCloud3.Location(Keypoint_mesh3.LocationIndex,:);
Keypoint_mesh3.Scale = ScaleParameter3D3(Keypoint_mesh3.Level);
Keypoint_mesh3.Count = length(Keypoint_mesh3.LocationIndex)


LevelTable = tabulate(Keypoint_mesh3.Level);
figure
plot(LevelTable(:,1)./PointCloud3.Resolution, ...
     LevelTable(:,2)./max(LevelTable(:,2)), 'k-', 'Linewidth', 4)


ZeroLogic3 = Keypoint_mesh3.Scale < 3*PointCloud3.Resolution;
         
FNames = fieldnames(Keypoint_mesh3);
for jj = 1 : length(FNames)
    if strcmpi(FNames{jj},'Count')
        Keypoint_mesh3 = rmfield(Keypoint_mesh3, FNames{jj});
    else
        Keypoint_mesh3.(FNames{jj})(ZeroLogic3,:) = [];
    end
end
Keypoint_mesh3.Count = length(Keypoint_mesh3.LocationIndex)



% NMSKeypoint_mesh = applyNMS(PointCloud3, DoG_mesh3, Keypoint_mesh3, t_scale, t_range1, 'sigma', DoGNormalize3D, CompareMethod);


% figure
% for i = 1 : NumSteps    
%     a = reshape(Signal_mesh(:,i), ImageSize(1), ImageSize(2));
%     imshow(a)    
% end


[SortValue3D, SortOrder3D] = sort(Keypoint_mesh3.Scale,'descend');


CircleAngle = linspace(0, 2*pi, 360);
CircleX = cos(CircleAngle);
CircleY = sin(CircleAngle);


PlotFeatures = 1:Keypoint_mesh3.Count;


figure
imshow(reshape(PointCloud.Signal, ImageSize(1), ImageSize(2)))
hold on
title(sprintf('Mesh: t=%0.2f',tau3))

for i = 1 : length(PlotFeatures)
% % % % %     CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY], [PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),1);PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),2)]);
    
    CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY]/beta,...
        [Keypoint_mesh3.Location(SortOrder3D(PlotFeatures(i)),2)/beta;...
        Keypoint_mesh3.Location(SortOrder3D(PlotFeatures(i)),1)/beta]);


    plot(CircleAtPoint(2,:), CircleAtPoint(1,:),'r')
    

%     pause
end









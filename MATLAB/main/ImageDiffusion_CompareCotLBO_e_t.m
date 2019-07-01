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
    




setTau = @(e) e;

tau = setTau(PointCloud.Resolution);

ScaleParameter3D = findScaleParameter(tau, 1, NumSteps, 'Laplacian', 'Natural');


[ItL, LBM] = makeCotangentLaplaceBeltrami( ImageFileNameOff, BDF, tau, alpha);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Diffusion w/ umb LBO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);

DoG = buildDoG(Signal, ScaleParameter3D, 'NLoG'); % DoGNormalize3D


[Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors3D.Distance);
Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
Keypoint.Scale = ScaleParameter3D(Keypoint.Level);
Keypoint.Count = length(Keypoint.LocationIndex)


LevelTable = tabulate(Keypoint.Level);
figure
plot(LevelTable(:,1)./PointCloud.Resolution, ...
     LevelTable(:,2)./max(LevelTable(:,2)), 'k-', 'Linewidth', 4)


ZeroLogic = Keypoint.Scale < 3*PointCloud.Resolution;
         
FNames = fieldnames(Keypoint);
for jj = 1 : length(FNames)
    if strcmpi(FNames{jj},'Count')
        Keypoint = rmfield(Keypoint, FNames{jj});
    else
        Keypoint.(FNames{jj})(ZeroLogic,:) = [];
    end
end
Keypoint.Count = length(Keypoint.LocationIndex)



% NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range1, 'sigma', DoGNormalize3D, CompareMethod);
% figure
% for i = 1 : NumSteps    
%     a = reshape(Signal(:,i), ImageSize(1), ImageSize(2));
%     imshow(a)    
% end


[SortValue3D, SortOrder3D] = sort(Keypoint.Scale,'descend');


CircleAngle = linspace(0, 2*pi, 360);
CircleX = cos(CircleAngle);
CircleY = sin(CircleAngle);


PlotFeatures = 1:Keypoint.Count;%100000:101000;%Keypoint.Count;


figure
imshow(reshape(PointCloud.Signal, ImageSize(1), ImageSize(2)))
hold on
title(sprintf('Mesh: t=%0.2f',tau))

% Still need to work on. 1/3/2019
for i = 1 : length(PlotFeatures)
% % % % %     CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY], [PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),1);PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),2)]);
    
    CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY],...
        [Keypoint.Location(SortOrder3D(PlotFeatures(i)),2);...
        Keypoint.Location(SortOrder3D(PlotFeatures(i)),1)]);


    plot(CircleAtPoint(2,:), CircleAtPoint(1,:),'r')
    

%     pause
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta = 3;

PointCloud3 = PointCloud;
PointCloud3.Location = PointCloud3.Location * beta;
PointCloud3 = findLocalResolution(PointCloud3, Neighbors3D.Distance);
PointCloud3 = findMeshResolution(PointCloud3, 'Model')

tau3 = setTau(PointCloud3.Resolution);

[ItL3, LBM3] = makeCotangentLaplaceBeltrami( ImageFileNameOff3, BDF, tau3, alpha);



% ItL3{1,1} = speye(PointCloud.LocationCount,PointCloud.LocationCount) - beta^2*tau*LBM3;



% ScaleParameter3D3 = findScaleParameter(beta^2*tau, alpha, beta*NumSteps, 'Laplacian', 'Natural');
ScaleParameter3D3 = findScaleParameter(tau3, alpha, beta*NumSteps, 'Laplacian', 'Natural');


Signal3 = performBDFDiffusion(PointCloud3.Signal, beta*NumSteps, ItL3);


DoG3 = buildDoG(Signal3, ScaleParameter3D3, 'NLoG'); % DoGNormalize3D


[Keypoint3.LocationIndex, Keypoint3.Level] = findKeypoint_c(DoG3, Neighbors3D.Distance);
Keypoint3.Location = PointCloud3.Location(Keypoint3.LocationIndex,:);
Keypoint3.Scale = ScaleParameter3D3(Keypoint3.Level);
Keypoint3.Count = length(Keypoint3.LocationIndex)


LevelTable = tabulate(Keypoint3.Level);
figure
plot(LevelTable(:,1)./PointCloud3.Resolution, ...
     LevelTable(:,2)./max(LevelTable(:,2)), 'k-', 'Linewidth', 4)


ZeroLogic3 = Keypoint3.Scale < 3*PointCloud3.Resolution;
         
FNames = fieldnames(Keypoint3);
for jj = 1 : length(FNames)
    if strcmpi(FNames{jj},'Count')
        Keypoint3 = rmfield(Keypoint3, FNames{jj});
    else
        Keypoint3.(FNames{jj})(ZeroLogic3,:) = [];
    end
end
Keypoint3.Count = length(Keypoint3.LocationIndex)



% NMSKeypoint = applyNMS(PointCloud3, DoG3, Keypoint3, t_scale, t_range1, 'sigma', DoGNormalize3D, CompareMethod);


% figure
% for i = 1 : NumSteps    
%     a = reshape(Signal(:,i), ImageSize(1), ImageSize(2));
%     imshow(a)    
% end


[SortValue3D, SortOrder3D] = sort(Keypoint3.Scale,'descend');


CircleAngle = linspace(0, 2*pi, 360);
CircleX = cos(CircleAngle);
CircleY = sin(CircleAngle);


PlotFeatures = 1:Keypoint3.Count;


figure
imshow(reshape(PointCloud.Signal, ImageSize(1), ImageSize(2)))
hold on
title(sprintf('Mesh: t=%0.2f',tau3))

for i = 1 : length(PlotFeatures)
% % % % %     CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY], [PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),1);PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),2)]);
    
    CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY]/beta,...
        [Keypoint3.Location(SortOrder3D(PlotFeatures(i)),2)/beta;...
        Keypoint3.Location(SortOrder3D(PlotFeatures(i)),1)/beta]);


    plot(CircleAtPoint(2,:), CircleAtPoint(1,:),'r')
    

%     pause
end











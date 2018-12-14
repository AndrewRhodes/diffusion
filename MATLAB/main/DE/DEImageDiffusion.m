

close all
clear
clc

global ProjectRoot; % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Model = 'Sunflower_e1_50445';


alpha = 1;
BDF = 1;
NumSteps = 2000;

options.rho = 6;
options.dtype = 'euclidean'; % 'euclidean', 'geodesic' %
options.htype = 'ddr'; % 'psp', 'ddr'

DoGNormalize = 'NLoG'; % 'DoG', 'AbsDoG', 'NLoG', 'AbsNLoG'
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'


setTau = @(e_bar) (e_bar/8)^2;
setHs = @(e_bar) 2*e_bar^(1/5);


t_scale = 0.7;
t_range = 1/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FileNameImage = strcat(ProjectRoot,'/models/sunflowers_285x177.jpeg');
% FileNameImage = strcat(ProjectRoot,'/models/Comet67P407x564.png');

Image = double(rgb2gray(imread(FileNameImage)))./255;

[ImageHeight, ImageWidth] = size(Image);
NumPixels = ImageHeight * ImageWidth;


FileNameOff = strcat(ProjectRoot,'/models/',Model,'.off');

if ~exist(FileNameOff, 'file')
    [xSurf3D, ySurf3D, zSurf3D] = ndgrid(1:ImageHeight,1:ImageWidth,0);
    PointCloud.Location = [ySurf3D(:), xSurf3D(:), zSurf3D(:)];
    PointCloud.Face = delaunay(xSurf3D(:), ySurf3D(:));
    
    save_off(PointCloud.Location, PointCloud.Face, FileNameOff);
else
    [PointCloud.Location, PointCloud.Face] = read_off(FileNameOff);
    
    PointCloud.Location = PointCloud.Location';
    PointCloud.Face = PointCloud.Face';    
end

PointCloud.LocationCount = length(PointCloud.Location);
PointCloud.FaceCount = length(PointCloud.Face);
PointCloud.Signal = Image(:);
PointCloud = findMeshResolution(PointCloud, 'Model');



FileNameNeighbors = strcat(ProjectRoot, '/models/Neighbors_',Model,'.mat');

if ~exist( FileNameNeighbors, 'file')
    [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);
    save(FileNameNeighbors ,'Neighbors', '-v7.3')
else
    load(FileNameNeighbors, 'Neighbors');
end

PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau = setTau(PointCloud.Resolution);

if strcmp(options.htype, 'psp')
    options.hs = setHs(PointCloud.Resolution);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileNameItL = strcat(ProjectRoot,'/models/',Model,'_ItL','_BDF',num2str(BDF),...
                '_a',num2str(alpha),'_tau',num2str(tau),...
                '_rho',num2str(options.rho),'_',options.dtype(1:3),'_',...
                options.htype,'.mat');


% if ~exist( FileNameItL,'file')
    ItL = makeMeshLaplaceBeltrami( FileNameOff, options, BDF, tau, alpha);
%     save(FileNameItL, 'ItL', '-v7.3');
% else
%     load(FileNameItL, 'ItL');
% end
        

ScaleParameter = findScaleParameter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);
        
DoG = buildDoG(Signal, ScaleParameter, DoGNormalize);


[Keypoint.LocationIndex, Keypoint.Level] = findKeypoint_c(DoG, Neighbors.Distance);
% Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
Keypoint.Scale = ScaleParameter(Keypoint.Level);
Keypoint.Count = length(Keypoint.LocationIndex);

NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod);
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Images to test scale and location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = linspace(0, 2*pi, 300);
CircleX = cos(theta');
CircleY = sin(theta');
% 
% [SphereX, SphereY, SphereZ] = sphere(100);
% SphereX = reshape(SphereX,[],1);
% SphereY = reshape(SphereY,[],1);
% SphereZ = reshape(SphereZ,[],1);


[Value, LocationIndex] = sort(Keypoint.Scale,'descend');

figure
imshow(Image);
hold on

for i = 1 : 2000 %Keypoint.Count
    CircleAtPoint = bsxfun(@plus, Keypoint.Scale(LocationIndex(i))*[CircleX, CircleY], ...
        Keypoint.Location(LocationIndex(i),1:2));

    plot(CircleAtPoint(:,1),CircleAtPoint(:,2),'r.')
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma2D = 0.75;

% ScaleParameter2D = findScaleParameter(sqrt(2*tau), 1, NumSteps, 'Gaussian', 'Natural');
% Gsigma = sqrt(2*tau);
ScaleParameter2D = findScaleParameter(sigma2D, 1, NumSteps, 'Gaussian', 'Natural');
Gsigma = sigma2D;

gausswindow = 5;
% [x,y] = meshgrid(-gausswindow:gausswindow,-gausswindow:gausswindow);
% G =  (1/(2*pi*Gsigma^2)) .* exp( -(x.^2 + y.^2) ./ (2*Gsigma^2) );
% G = G./sum(G(:));

Gauss = fspecial('gaussian', [gausswindow, gausswindow], Gsigma);

Signal2D = zeros(ImageHeight, ImageWidth, NumSteps);
Signal2D(:,:,1) = Image;

% figure
for i = 1 : NumSteps - 1
   Signal2D(:,:,i+1) = imfilter(Signal2D(:,:,i), Gauss, 'replicate', 'same', 'conv');
%    Signal2D(:,:,i+1) = Signal2D(:,:,i+1) ./ max(max(Signal2D(:,:,i+1)));
%    imshow(Signal2D(:,:,i+1))
end


DoG2D = zeros(ImageHeight, ImageWidth, NumSteps-1);

% figure
for i = 1 : NumSteps - 1
   DoG2D(:,:,i) = Signal2D(:,:,i+1) - Signal2D(:,:,i); 
   DoG2D(:,:,i) = DoG2D(:,:,i) ./ max(max(DoG2D(:,:,i)));
%    imshow(DoG2D(:,:,i))
end
DoG2DReshape = reshape(DoG2D, NumPixels, NumSteps-1);

Neighbors2D = cell(NumPixels,1);
for i = 1 : NumPixels
    
    [CurX, CurY] = ind2sub([ImageHeight,ImageWidth], i);
    CurNeigh = [CurX-1, CurY-1;
                CurX-1, CurY;
                CurX-1, CurY+1;
                CurX+1, CurY+1;
                CurX+1, CurY;
                CurX+1, CurY-1;
                CurX, CurY-1;
                CurX, CurY+1;];
%     AboveNeig = [CurNeigh; CurX, CurY];
%     BelowNeigh = [CurNeigh; CurX, CurY];
    
    if ~any(CurNeigh(:) == 0) && ~any(CurNeigh(:,1) > ImageHeight) && ~any(CurNeigh(:,2) > ImageWidth)
        
        CurNeigh = sub2ind([ImageHeight,ImageWidth], CurNeigh(:,1), CurNeigh(:,2));
        
        Neighbors2D{i,1} = [CurNeigh];
    else
        Neighbors2D{i,1} = 0;
    end
    
end

[Keypoint2D.LocationIndex, Keypoint2D.Level] = findKeypoint_c(DoG2DReshape, Neighbors2D);
Keypoint2D.Scale = ScaleParameter2D(Keypoint2D.Level);
[a,b] = ind2sub([ImageHeight,ImageWidth], Keypoint2D.LocationIndex);
Keypoint2D.Location = [a,b];
Keypoint2D.Count = length(Keypoint2D.LocationIndex);




% RemoveTooSmall = (Keypoint2D.Level < 10);
% 
% if isfield(Keypoint2D, 'Count')
%     Keypoint2D = rmfield(Keypoint2D, 'Count');
% end
% FNames = fieldnames(Keypoint2D);
% for i = 1 : length(FNames)
%     Keypoint2D.(FNames{i})(RemoveTooSmall,:) = [];
% end
% Keypoint2D.Count = length(Keypoint2D.LocationIndex);





theta = linspace(0, 2*pi, 300);
CircleX = cos(theta');
CircleY = sin(theta');


[Value, LocationIndex] = sort(Keypoint2D.Scale,'descend');


figure
imshow(Image);
hold on

for i = 1 : 2000 % Keypoint2D.Count
    
    CircleAtPoint = bsxfun(@plus, Keypoint2D.Scale(LocationIndex(i))*[CircleX, CircleY], ...
        Keypoint2D.Location(LocationIndex(i),1:2));

    plot(CircleAtPoint(:,2),CircleAtPoint(:,1),'r.')
end

































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
% NumSteps = 500;
FileLocation = strcat(ProjectRoot,'/models/image/');
ImageChoice = 'Woodburn' %'Voorhees' % 'Woodburn' % 'sunflowers' % 'comet67p' % 



% Settings: 2D image diffusion
gausswindow = 5;
sigma2D = sqrt(2); % 0.75; % 
DoGNormalize2D = 'DoG';



% Settings: 3D surface diffusion
DoGNormalize3D = 'DoG';
CompareMethod = '<>'; % '<', '>', '<>'
KeypointMethod = 'Old'; % 'Old', 'New'
alpha = 1;
tau_n = sigma2D^2 / (2*alpha);
BDF = 1;

% l_scale = 1/sqrt(2);
l_range1 = 1/2;
l_range2 = 2;



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
   
elseif strcmpi(ImageChoice, 'Voorhees')
   
    ImageName = 'VoorheesComputingCenter_831x436.jpg';
    Image = double(rgb2gray(imread(strcat(FileLocation,ImageName))));
    Image = imresize(Image,0.5);
    Image = Image./255;
    
elseif strcmpi(ImageChoice, 'Woodburn')
    ImageName = 'WVUWoodburn_1095x709.jpg';
    ImageName = 'WVUWoodburn_466x303.jpg'; % 'WVUWoodburn_524x358.jpg';
    Image = double(rgb2gray(imread(strcat(FileLocation,ImageName))));
    Image = Image./255;
    
elseif strcmpi(ImageChoice, 'Woodburn2')

    ImageName = 'WVUWoodburn_438x284.jpg';
    Image = double(rgb2gray(imread(strcat(FileLocation,ImageName))));
    Image = Image./255;    

end


figure
imshow(Image)


ImageSize = size(Image);
NumVertex = numel(Image);
NumPixels = numel(Image);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



ImageFileName = strcat(ImageChoice ,'_e1_', num2str(NumPixels));
ImageFileName3 = strcat(ImageChoice ,'_e3_', num2str(NumPixels));
ImageFileNameOff = strcat(FileLocation, ImageFileName, '.off');
ImageFileNameOff3 = strcat(FileLocation, ImageFileName, '.off');

FileNameItL_mesh = strcat(ImageFileName,'_ItL','_BDF',num2str(BDF),'_rho',...
              num2str(options.rho),'_',options.dtype(1:3),'_',...
              options.htype,'_t0.25','_a',num2str(alpha),'.mat'); 

FileNameLBM_mesh = strcat(ImageFileName,'_LBM','_rho',num2str(options.rho),...
                   '_',options.dtype(1:3),'_',options.htype,'.mat'); 
               

FileNameItL_cot = strcat(ImageFileName,'_ItL','_BDF',num2str(BDF),'_cot'...
                  ,'_t0.25','_a',num2str(alpha),'.mat'); 

FileNameLBM_cot = strcat(ImageFileName,'_LBM','_cot','.mat'); 
               



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Image Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% c = 0;
% N = 9
% tn = (sqrt(2)).^(2*(0:N-1)')
% tau_n = [diff(tn)]
% ScaleParameter2D =  sqrt( 2 * alpha * tn)
% sigma_n = sqrt(diff(ScaleParameter2D.^2))
% 
% Signal2D = zeros(ImageSize(1), ImageSize(2), N);
% Signal2D(:, :, 1) = Image;
% 
% % ScaleParamterImage1 = zeros(NumSteps,1);
% 
% figure
% for i = 1 : N - 1
%     gausswindow = max(5, ceil(3*sigma_n(i)));
%     gausswindow = gausswindow + mod(gausswindow+1,2)
%     
%     Gauss = fspecial('gaussian', [gausswindow, gausswindow], sigma_n(i));
%     
%     Signal2D(:,:,i+1) = imfilter(Signal2D(:,:,i), Gauss, 'replicate');
% 
%     imshow(Signal2D(:,:,i+1))
%     drawnow
% %     pause
%     
% %     ScaleParameter2D(i+1) = sqrt(i) * sigma2D;
%     
% end
% 
% DiffusedSignal2D = reshape(Signal2D, NumPixels, N);
% 
% Neighbors2D = findAdjacentNeighbors2D(ImageSize);
% 
% DoG2D = buildDoG(DiffusedSignal2D, ScaleParameter2D, 'DoG');
% 
% % for i = 1 : NumSteps - 1
% %     imshow(reshape(DoG2D(:,i), ImageSize(1), ImageSize(2)),[])
% %     pause
% % end
% 
% [Keypoint2D.LocationIndex, Keypoint2D.Level] = findKeypoint_c(DoG2D, Neighbors2D);
% Keypoint2D.Scale = ScaleParameter2D(Keypoint2D.Level);
% [Iindex, Jindex] = ind2sub([ImageSize(1),ImageSize(2)], Keypoint2D.LocationIndex);
% Keypoint2D.Location = [Iindex, Jindex];
% Keypoint2D.Count = length(Keypoint2D.LocationIndex)
% 
% 
% Remove2D = false(size(Keypoint2D.Scale));
% for i = 1 : Keypoint2D.Count
%     if abs(DoG2D(Keypoint2D.LocationIndex(i), Keypoint2D.Level(i))) < 0.03
%         Remove2D(i) = 1;
%     end
% end
% 
% 
% ZeroLogic = Keypoint2D.Scale < 3 ;
%          
% f = 1.5;
% Remove = false(size(Keypoint2D.Scale));
% for i = 1 : Keypoint2D.Count
%     
%     if (~Remove(i)) && (Keypoint2D.Location(i,1) < f*Keypoint2D.Scale(i))
%         Remove(i) = 1;
%     end
%     
%     if (~Remove(i)) && (Keypoint2D.Location(i,1) > ImageSize(2)-f*Keypoint2D.Scale(i))
%         Remove(i) = 1;
%     end
%     
%     if (~Remove(i)) && (Keypoint2D.Location(i,2) < f*Keypoint2D.Scale(i))
%         Remove(i) = 1;
%     end
%     
%     if (~Remove(i)) && (Keypoint2D.Location(i,2) > ImageSize(1)-f*Keypoint2D.Scale(i))
%         Remove(i) = 1;
%     end
%     
% end
% nnz(Remove)
% 
% 
% FNames = fieldnames(Keypoint2D);
% for jj = 1 : length(FNames)
%     if strcmpi(FNames{jj},'Count')
%         Keypoint2D = rmfield(Keypoint2D, FNames{jj});
%     else
%         Keypoint2D.(FNames{jj})(Remove2D | Remove | ZeroLogic,:) = [];
%     end
% end
% Keypoint2D.Count = length(Keypoint2D.LocationIndex)
% 
% 
% 
% 
% LevelTable = tabulate(Keypoint2D.Level);
% figure
% plot(ScaleParameter2D(LevelTable(:,1)), ...
%      LevelTable(:,2)./max(LevelTable(:,2)), 'k-', 'Linewidth', 4)
% 
%  
% NMSKeypoint2D = applyNMS(PointCloud, DoG2D, Keypoint2D, 1/sqrt(2), l_range1, 'sigma', 'DoG', CompareMethod)
% % NMSKeypoint2D = applyNMS(PointCloud, DoG2D, Keypoint2D, 1/sqrt(2), l_range2, 'ebar', 'DoG', CompareMethod)
% 
% 
% [SortValue2D, SortOrder2D] = sort(NMSKeypoint2D.Scale,'ascend');
% 
% CircleAngle = linspace(0, 2*pi, 360);
% CircleX = cos(CircleAngle);
% CircleY = sin(CircleAngle);
% 
% PlotFeatures = 1:1000;%5:max(SortOrder2D);
% 
% figure
% imshow(Image)
% hold on
% 
% for i = 1 : length(PlotFeatures)
% 
%     [a,b] = ind2sub([ImageSize(1),ImageSize(2)], NMSKeypoint2D.LocationIndex(SortOrder2D(PlotFeatures(i))));
%     
%     CircleAtPoint = bsxfun(@plus, SortValue2D(PlotFeatures(i))*[CircleX; CircleY], ...
%         [NMSKeypoint2D.Location(SortOrder2D(PlotFeatures(i)),1); NMSKeypoint2D.Location(SortOrder2D(PlotFeatures(i)),2)]);
%     
%     plot(CircleAtPoint(2,:), CircleAtPoint(1,:),'r')
%     
% %         drawnow
% 
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Image Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Must add SIFT path to MATLAB 
run('/media/andrew/OSDisk/vlfeat/toolbox/vl_setup.m')
% Make sure it worked
vl_version verbose



[Keypoint_SIFT, Descriptor1] = vl_sift(single(Image), 'WindowSize', 3);
size(Keypoint_SIFT,2)


NMSKeypoint_SIFT = Keypoint_SIFT;
NMSKeypoint_SIFT(:, NMSKeypoint_SIFT(3,:)< 3) = [];
size(NMSKeypoint_SIFT,2)


CircleAngle = linspace(0, 2*pi, 360);
CircleX = cos(CircleAngle);
CircleY = sin(CircleAngle);


figure
imshow(Image)
hold on

for i = 1 : size(NMSKeypoint_SIFT,2)

%     [a,b] = ind2sub([ImageSize(1),ImageSize(2)], Keypoint2D.LocationIndex(SortOrder2D(PlotFeatures(i))));
    
    CircleAtPoint = bsxfun(@plus, NMSKeypoint_SIFT(3,i)*[CircleX; CircleY], ...
        [NMSKeypoint_SIFT(1,i); NMSKeypoint_SIFT(2,i)]);
    
    plot(CircleAtPoint(1,:), CircleAtPoint(2,:),'r', 'linewidth', 3)
    
%         drawnow

end




% Signal2D = zeros(ImageSize(1), ImageSize(2), NumSteps);
% Signal2D(:, :, 1) = Image;
% 
% Gauss = fspecial('gaussian', [gausswindow, gausswindow], sigma2D);
% ScaleParameter2D = zeros(NumSteps,1);
% ScaleParamterImage1 = zeros(NumSteps,1);
% 
% % figure
% for i = 1 : NumSteps - 1
%         
%     Signal2D(:,:,i+1) = imfilter(Signal2D(:,:,i), Gauss, 'replicate');
% 
% %     imshow(Signal2D(:,:,i+1))
% %     drawnow
% %     pause
%     
%     ScaleParameter2D(i+1) = sqrt(i) * sigma2D;
%     
% end
% 
% DiffusedSignal2D = reshape(Signal2D, NumPixels, NumSteps);
% 
% Neighbors2D = findAdjacentNeighbors2D(ImageSize);
% 
% DoG2D = buildDoG(DiffusedSignal2D, ScaleParameter2D, 'DoG');
% 
% % for i = 1 : NumSteps - 1
% %     imshow(reshape(DoG2D(:,i), ImageSize(1), ImageSize(2)),[])
% %     pause
% % end
% 
% [Keypoint2D.LocationIndex, Keypoint2D.Level] = findKeypoint_c(DoG2D, Neighbors2D);
% Keypoint2D.Scale = ScaleParameter2D(Keypoint2D.Level);
% [Iindex, Jindex] = ind2sub([ImageSize(1),ImageSize(2)], Keypoint2D.LocationIndex);
% Keypoint2D.Location = [Iindex, Jindex];
% Keypoint2D.Count = length(Keypoint2D.LocationIndex)
% 
% LevelTable = tabulate(Keypoint2D.Level);
% figure
% plot(ScaleParameter2D(LevelTable(:,1)), ...
%      LevelTable(:,2)./max(LevelTable(:,2)), 'k-', 'Linewidth', 4)
% xlim([0 5])
% 
% 
% % Comet67P Only
% ZeroLogic = (Image(Keypoint2D.LocationIndex) < 0.04) | ...
%              (Keypoint2D.Location(:,2) <= 3) | ...
%              (Keypoint2D.Location(:,1) <= 3) | ...
%              (Keypoint2D.Location(:,1) >= ImageSize(1)-3) | ...
%              (Keypoint2D.Location(:,2) >= ImageSize(2)-3) | ...
%              (Keypoint2D.Scale < 3);
%          
% FNames = fieldnames(Keypoint2D);
% for jj = 1 : length(FNames)
%     if strcmpi(FNames{jj},'Count')
%         Keypoint2D = rmfield(Keypoint2D, FNames{jj});
%     else
%         Keypoint2D.(FNames{jj})(ZeroLogic,:) = [];
%     end       
% end
% Keypoint2D.Count = length(Keypoint2D.LocationIndex)
% 
% 
% 
% [SortValue2D, SortOrder2D] = sort(Keypoint2D.Scale,'ascend');
% 
% CircleAngle = linspace(0, 2*pi, 360);
% CircleX = cos(CircleAngle);
% CircleY = sin(CircleAngle);
% 
% 
% PlotFeatures = 10000:11000;%5:max(SortOrder2D);
% 
% figure
% imshow(Image)
% hold on
% 
% for i = 1 : length(PlotFeatures)
% 
%     [a,b] = ind2sub([ImageSize(1),ImageSize(2)], Keypoint2D.LocationIndex(SortOrder2D(PlotFeatures(i))));
%     
%     CircleAtPoint = bsxfun(@plus, SortValue2D(PlotFeatures(i))*[CircleX; CircleY], ...
%         [Keypoint2D.Location(SortOrder2D(PlotFeatures(i)),1); Keypoint2D.Location(SortOrder2D(PlotFeatures(i)),2)]);
%     
%     plot(CircleAtPoint(2,:), CircleAtPoint(1,:),'r')
%     
% %         drawnow
% 
% end




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

         
% if ~exist(strcat(FileLocation, FileNameItL_mesh), 'file')
    [LBM_mesh] = makeMeshLaplaceBeltrami( ImageFileNameOff, options);
%     [ItL, LBM, A, W] = makeMeshLaplaceBeltrami( ImageFileNameOff, options, BDF, tau_n, alpha);
%     ItL_mesh = ItL;
%     LBM_mesh = LBM;
%     save(strcat(FileLocation, FileNameItL_mesh),'ItL','-v7.3');
%     save(strcat(FileLocation, FileNameLBM_mesh),'LBM','-v7.3');
%     clear ItL LBM
% else
%     load(strcat(FileLocation, FileNameItL_mesh), 'ItL');
%     load(strcat(FileLocation, FileNameLBM_mesh), 'LBM');
%     ItL_mesh = ItL;
%     LBM_mesh = LBM;
%     clear ItL LBM
% end


% Signal_mesh = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL_mesh);


l_ebar = ceil(max(Keypoint_SIFT(3,:)));
k = sqrt(2);
l_scale = 1/k;
set_t0 = @(e_bar) e_bar ;
set_ltau = @(e_bar) 2*e_bar;
set_NumSteps = @(e_bar, t_0) ceil(log((l_ebar*e_bar)^2 / (2*alpha*t_0)) / (2*log(k)));


l_tau = set_ltau(PointCloud.Resolution);
t_0 = set_t0(PointCloud.Resolution);
N = set_NumSteps(PointCloud.Resolution, t_0)

[tn, tau_n, sigma_n] = findScaleStep(k, t_0, alpha, N);

[Signal_mesh, IterCount] = performBDFDiffusion(PointCloud.Signal, LBM_mesh, alpha, tau_n, N, l_tau);
    

DoG_mesh = buildDoG(Signal_mesh, sigma_n, DoGNormalize3D); % DoGNormalize3D


[Keypoint_mesh.LocationIndex, Keypoint_mesh.Level] = findKeypoint_c(DoG_mesh, Neighbors3D.Distance);
Keypoint_mesh.Location = PointCloud.Location(Keypoint_mesh.LocationIndex,:);
Keypoint_mesh.Scale = sigma_n(Keypoint_mesh.Level);
Keypoint_mesh.Count = length(Keypoint_mesh.LocationIndex)


LevelTable = tabulate(Keypoint_mesh.Level);
figure
plot(sigma_n(LevelTable(:,1))/PointCloud.Resolution, ...
     LevelTable(:,2)./max(LevelTable(:,2)), 'k-', 'Linewidth', 4)

 
 
NMSKeypoint_mesh = applyNMS(PointCloud, DoG_mesh, Keypoint_mesh, l_scale, l_range1, 'sigma', DoGNormalize3D, CompareMethod)
% NMSKeypoint_mesh = applyNMS(PointCloud, DoG_mesh, Keypoint_mesh, l_scale, l_range2, 'ebar', DoGNormalize3D, CompareMethod)


 

ZeroLogic = NMSKeypoint_mesh.Scale < 3 * PointCloud.Resolution;
         

f = 1.5;
Remove = false(size(NMSKeypoint_mesh.Scale));
for i = 1 : NMSKeypoint_mesh.Count
    
    if (~Remove(i)) && (NMSKeypoint_mesh.Location(i,1) < f*NMSKeypoint_mesh.Scale(i))
        Remove(i) = 1;
    end
    
    if (~Remove(i)) && (NMSKeypoint_mesh.Location(i,1) > ImageSize(2)-f*NMSKeypoint_mesh.Scale(i))
        Remove(i) = 1;
    end
    
    if (~Remove(i)) && (NMSKeypoint_mesh.Location(i,2) < f*NMSKeypoint_mesh.Scale(i))
        Remove(i) = 1;
    end
    
    if (~Remove(i)) && (NMSKeypoint_mesh.Location(i,2) > ImageSize(1)-f*NMSKeypoint_mesh.Scale(i))
        Remove(i) = 1;
    end
    
end
nnz(Remove)


FNames = fieldnames(NMSKeypoint_mesh);
for jj = 1 : length(FNames)
    if strcmpi(FNames{jj},'Count')
        NMSKeypoint_mesh = rmfield(NMSKeypoint_mesh, FNames{jj});
    else
        NMSKeypoint_mesh.(FNames{jj})(ZeroLogic | Remove,:) = [];
    end
end
NMSKeypoint_mesh.Count = length(NMSKeypoint_mesh.LocationIndex)


% LevelTable = tabulate(Keypoint_mesh.Level);
% figure
% plot(LevelTable(:,1), ...
%      LevelTable(:,2)./max(LevelTable(:,2)), 'k-', 'Linewidth', 4)
  
 

% figure
% for i = 1 : NumSteps    
%     a = reshape(Signal_mesh(:,i), ImageSize(1), ImageSize(2));
%     imshow(a)    
% end


[SortValue3D, SortOrder3D] = sort(NMSKeypoint_mesh.Scale,'descend');


CircleAngle = linspace(0, 2*pi, 360);
CircleX = cos(CircleAngle);
CircleY = sin(CircleAngle);


PlotFeatures = 1:NMSKeypoint_mesh.Count;%:max(SortOrder3D);


figure
imshow(reshape(PointCloud.Signal, ImageSize(1), ImageSize(2)))
hold on
% title('Mesh')

for i = 1 : length(PlotFeatures)
% % % % %     CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY], [PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),1);PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),2)]);
    
    CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY],...
        [NMSKeypoint_mesh.Location(SortOrder3D(PlotFeatures(i)),2);...
        NMSKeypoint_mesh.Location(SortOrder3D(PlotFeatures(i)),1)]);


    plot(CircleAtPoint(2,:), CircleAtPoint(1,:),'r', 'linewidth', 3)
    

%     pause
end









%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Diffusion w/ cotangent LBO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% if ~exist( strcat(FileLocation, FileNameItL_cot), 'file')
    LBM_cot = makeCotangentLaplaceBeltrami( ImageFileNameOff );
%     ItL_cot = ItL;
%     LBM_cot = LBM;
%     save( strcat(FileLocation, FileNameItL_cot), 'ItL', '-v7.3');
%     save( strcat(FileLocation, FileNameLBM_cot), 'LBM', '-v7.3');
%     clear ItL LBM
% else
%     load( strcat(FileLocation, FileNameItL_cot), 'ItL');
%     load( strcat(FileLocation, FileNameLBM_cot), 'LBM');
%     ItL_cot = ItL;
%     LBM_cot = LBM;
%     clear ItL LBM
% end



[Signal_cot, IterCount] = performBDFDiffusion(PointCloud.Signal, LBM_cot, alpha, tau_n, N, l_tau);
    

DoG_cot = buildDoG(Signal_cot, sigma_n, DoGNormalize3D);


[Keypoint_cot.LocationIndex, Keypoint_cot.Level] = findKeypoint_c(DoG_cot, Neighbors3D.Connect);
Keypoint_cot.Location = PointCloud.Location(Keypoint_cot.LocationIndex,:);
Keypoint_cot.Scale = sigma_n(Keypoint_cot.Level);
Keypoint_cot.Count = length(Keypoint_cot.Scale)



LevelTable = tabulate(Keypoint_cot.Level);
figure
plot(sigma_n(LevelTable(:,1))/sqrt(PointCloud.Resolution), ...
     LevelTable(:,2)./max(LevelTable(:,2)), 'k-', 'Linewidth', 4)

 
 
NMSKeypoint_cot = applyNMS(PointCloud, DoG_mesh, Keypoint_cot, l_scale, l_range1, 'sigma', DoGNormalize3D, CompareMethod)
% NMSKeypoint_cot = applyNMS(PointCloud, DoG_mesh, Keypoint_cot, l_scale, l_range2, 'ebar', DoGNormalize3D, CompareMethod)


 

ZeroLogic = NMSKeypoint_cot.Scale < 3 * PointCloud.Resolution;
         

f = 1.5;
Remove = false(size(NMSKeypoint_cot.Scale));
for i = 1 : NMSKeypoint_cot.Count
    
    if (~Remove(i)) && (NMSKeypoint_cot.Location(i,1) < f*NMSKeypoint_cot.Scale(i))
        Remove(i) = 1;
    end
    
    if (~Remove(i)) && (NMSKeypoint_cot.Location(i,1) > ImageSize(2)-f*NMSKeypoint_cot.Scale(i))
        Remove(i) = 1;
    end
    
    if (~Remove(i)) && (NMSKeypoint_cot.Location(i,2) < f*NMSKeypoint_cot.Scale(i))
        Remove(i) = 1;
    end
    
    if (~Remove(i)) && (NMSKeypoint_cot.Location(i,2) > ImageSize(1)-f*NMSKeypoint_cot.Scale(i))
        Remove(i) = 1;
    end
    
end
nnz(Remove)


FNames = fieldnames(NMSKeypoint_cot);
for jj = 1 : length(FNames)
    if strcmpi(FNames{jj},'Count')
        NMSKeypoint_cot = rmfield(NMSKeypoint_cot, FNames{jj});
    else
        NMSKeypoint_cot.(FNames{jj})(ZeroLogic | Remove,:) = [];
    end
end
NMSKeypoint_cot.Count = length(NMSKeypoint_cot.LocationIndex)




[SortValue3D, SortOrder3D] = sort(NMSKeypoint_cot.Scale,'descend');


CircleAngle = linspace(0, 2*pi, 360);
CircleX = cos(CircleAngle);
CircleY = sin(CircleAngle);


PlotFeatures = 1:NMSKeypoint_cot.Count;


figure
imshow(reshape(PointCloud.Signal, ImageSize(1), ImageSize(2)))
hold on
% title('Cotangent')

for i = 1 : length(PlotFeatures)
% % % % %     CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY], [PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),1);PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),2)]);
    
CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY],...
        [NMSKeypoint_cot.Location(SortOrder3D(PlotFeatures(i)),2);...
        NMSKeypoint_cot.Location(SortOrder3D(PlotFeatures(i)),1)]);

    plot(CircleAtPoint(2,:), CircleAtPoint(1,:),'r', 'linewidth', 3)
    

%     pause
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Diffusion w/ umbrella LBO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LBM_umb = makeUmbrellaLaplaceBeltrami(PointCloud, Neighbors3D.Connect);



[Signal_umb, IterCount] = performBDFDiffusion(PointCloud.Signal, LBM_umb, alpha, tau_n, N, l_tau);
    

DoG_umb = buildDoG(Signal_umb, sigma_n, DoGNormalize3D);


[Keypoint_umb.LocationIndex, Keypoint_umb.Level] = findKeypoint_c(DoG_umb, Neighbors3D.Connect);
Keypoint_umb.Location = PointCloud.Location(Keypoint_umb.LocationIndex,:);
Keypoint_umb.Scale = sigma_n(Keypoint_umb.Level)*PointCloud.Resolution;
Keypoint_umb.Count = length(Keypoint_umb.Scale)




NMSKeypoint_umb = applyNMS(PointCloud, DoG_umb, Keypoint_umb, l_scale, l_range1, 'sigma', DoGNormalize3D, CompareMethod)
% NMSKeypoint_umb = applyNMS(PointCloud, DoG_umb, Keypoint_umb, l_scale, l_range2, 'ebar', DoGNormalize3D, CompareMethod)



ZeroLogic = NMSKeypoint_umb.Scale < 3*PointCloud.Resolution;
         

f = 1;
Remove = false(size(NMSKeypoint_umb.Scale));
for i = 1 : NMSKeypoint_umb.Count
    
    if (~Remove(i)) && (NMSKeypoint_umb.Location(i,1) < f*NMSKeypoint_umb.Scale(i))
        Remove(i) = 1;
    end
    
    if (~Remove(i)) && (NMSKeypoint_umb.Location(i,1) > ImageSize(2)-f*NMSKeypoint_umb.Scale(i))
        Remove(i) = 1;
    end
    
    if (~Remove(i)) && (NMSKeypoint_umb.Location(i,2) < f*NMSKeypoint_umb.Scale(i))
        Remove(i) = 1;
    end
    
    if (~Remove(i)) && (NMSKeypoint_umb.Location(i,2) > ImageSize(1)-f*NMSKeypoint_umb.Scale(i))
        Remove(i) = 1;
    end
    
end
nnz(Remove)


FNames = fieldnames(NMSKeypoint_umb);
for jj = 1 : length(FNames)
    if strcmpi(FNames{jj},'Count')
        NMSKeypoint_umb = rmfield(NMSKeypoint_umb, FNames{jj});
    else
        NMSKeypoint_umb.(FNames{jj})(ZeroLogic | Remove,:) = [];
    end
end
NMSKeypoint_umb.Count = length(NMSKeypoint_umb.LocationIndex)




[SortValue3D, SortOrder3D] = sort(NMSKeypoint_umb.Scale,'descend');


CircleAngle = linspace(0, 2*pi, 360);
CircleX = cos(CircleAngle);
CircleY = sin(CircleAngle);


PlotFeatures = 1:NMSKeypoint_umb.Count;


figure
imshow(reshape(PointCloud.Signal, ImageSize(1), ImageSize(2)))
hold on
% title('Umbrella')


for i = 1 : length(PlotFeatures)
% % % % %     CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY], [PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),1);PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),2)]);
    CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY],...
        [NMSKeypoint_umb.Location(SortOrder3D(PlotFeatures(i)),2);...
        NMSKeypoint_umb.Location(SortOrder3D(PlotFeatures(i)),1)]);

    plot(CircleAtPoint(2,:), CircleAtPoint(1,:),'r', 'linewidth', 3)
    

%     pause
end





KeypointSaveLocation = strcat(ProjectRoot, '/data/imagekeypointdata/',ImageChoice,'/');

save( strcat(KeypointSaveLocation, 'Keypoint_mesh.mat'), 'Keypoint_mesh', '-v7.3') 
save( strcat(KeypointSaveLocation, 'Keypoint_cot.mat'), 'Keypoint_cot', '-v7.3') 
save( strcat(KeypointSaveLocation, 'Keypoint_umb.mat'), 'Keypoint_umb', '-v7.3') 
save( strcat(KeypointSaveLocation, 'Keypoint_SIFT.mat'), 'Keypoint_SIFT', '-v7.3') 

save( strcat(KeypointSaveLocation, 'NMSKeypoint_mesh.mat'), 'NMSKeypoint_mesh', '-v7.3') 
save( strcat(KeypointSaveLocation, 'NMSKeypoint_cot.mat'), 'NMSKeypoint_cot', '-v7.3') 
save( strcat(KeypointSaveLocation, 'NMSKeypoint_umb.mat'), 'NMSKeypoint_umb', '-v7.3') 
save( strcat(KeypointSaveLocation, 'NMSKeypoint_SIFT.mat'), 'NMSKeypoint_SIFT', '-v7.3') 














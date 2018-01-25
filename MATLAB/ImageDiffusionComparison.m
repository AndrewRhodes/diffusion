
% Andrew Rhodes
% ASEL
% 2017

% Perform scale space on an 2D image of an camera man and find key points.
% Sample the image onto a 3D explicit surface.
% Perform implicit surface diffusion, find key points, and compare.


close all
clear 
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional Paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd ~/Desktop/Ashish/'CS3 Code'/
addpath(genpath('~/Documents/Software/MeshLP/'))
addpath(genpath('~/Documents/Software/cp_matrices/'))
addpath('~/Desktop/Ashish/CS3 Code/')
addpath(genpath('~/GitProjects/pose/MATLAB_PointCloudDescriptors/OURCVFH/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ImageName = 'sunflowers.png';
% Image = double(rgb2gray(imread(ImageName)));
% Image = Image(:,200:500);
% Image = Image./255;



ImageName = 'Comet67PImage.png';
Image = double(imread(ImageName));
Image = imresize(Image,0.5);
Image = Image./255;




% imshow(Image)

% 2D image diffusion
% Max2DTau = 100;
Max2DTau = 600;
gausswindow = 5;
tau2D = 0.75;
MaxLevel = round(Max2DTau / tau2D);

MaxImageSize = size(Image);
MiddleImage = MaxImageSize / 2;
NumPixels = prod(MaxImageSize);


% 3D explcit surface diffusion
eSS = 1; % Explicit Surface Spacing
MaxSurfSize = size(Image);
NumVertex = prod(MaxSurfSize);


% 3D implicit surface
spacing = 0.5;
porder = 3;
Lorder = 2;
dim = 3;
tauImplicit = spacing / 4;
tauImplicit = tau2D^2/2;
MaxTauImplicit = 200;
NumStepsImplicit = round(MaxTauImplicit / tauImplicit);
NumStepsImplicit = MaxLevel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Impulse Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Signal2D = zeros(MaxImageSize(1), MaxImageSize(2), MaxLevel);
Signal2D(:, :, 1) = Image;

Gauss = fspecial('gaussian', [gausswindow, gausswindow], tau2D);
ScaleParamterImage = zeros(MaxLevel,1);
figure
for i = 1 : MaxLevel - 1
    
%     imshow(Signal2D(:,:,i),[])
%     pause
    
    Signal2D(:,:,i+1) = imfilter(Signal2D(:,:,i), Gauss,'replicate');

    imshow(Signal2D(:,:,i+1))
    drawnow
    i
%     ScaleParamterImage(i+1) = sqrt(2*i * tau2D^2);
    ScaleParamterImage(i+1) = sqrt(i) * tau2D;
end

% sigma^2 = [12.5 25 50]


for i = [10, 50, 100, 200, MaxLevel]
    
    
    figure
    imshow(Signal2D(:,:,i))
    
    
end



[Iindex, Jindex] = ind2sub([MaxImageSize(1),MaxImageSize(2)], (1:NumPixels)');

Neighbors2D = cell(NumPixels,1);

for i = 1 : NumPixels
    
    CurrentPoint = [Iindex(i), Jindex(i)];
    
    PossibleNeighbors = [Iindex(i)+1, Jindex(i);
                        Iindex(i)-1, Jindex(i);
                        Iindex(i), Jindex(i)+1;
                        Iindex(i), Jindex(i)-1;
                        Iindex(i)+1, Jindex(i)+1;
                        Iindex(i)+1, Jindex(i)-1;
                        Iindex(i)-1, Jindex(i)+1;
                        Iindex(i)-1, Jindex(i)-1];
                    
    PossibleNeighbors(PossibleNeighbors(:,1) > MaxImageSize(1), :) = [];
    PossibleNeighbors(PossibleNeighbors(:,1) < 1, :) = [];
    PossibleNeighbors(PossibleNeighbors(:,2) > MaxImageSize(2), :) = [];
    PossibleNeighbors(PossibleNeighbors(:,2) < 1, :) = [];
    
    Neighbors2D{i,1} = PossibleNeighbors;
    
end



% LaplacianScale = zeros(NumPixels, MaxLevel-1);
LaplacianScaleNormalized = zeros(NumPixels, MaxLevel-1);
% LaplaceScaleInvariant = zeros(NumPixels, MaxLevel-1);



for i = 1 : MaxLevel-1
    
    
    % Estimated Laplacian
%     LaplacianScale(:,i) = ( reshape(Signal2D(:,:,i+1),[],1) - reshape(Signal2D(:,:,i),[],1) ) / (ScaleParamterImage(i+1)^2-ScaleParamterImage(i)^2);
    
    % Scale Normalized Laplacian
    LaplacianScaleNormalized(:,i) = 2 * ( reshape(Signal2D(:,:,i+1),[],1) - reshape(Signal2D(:,:,i),[],1) ) * (ScaleParamterImage(i)^2 / (ScaleParamterImage(i+1)^2 - ScaleParamterImage(i)^2) );

    % Scale Invariant Laplace
%     Fbar = (sum(reshape(Signal2D(:,:,i+1),[],1)) / (NumPixels)) * ones(NumPixels, 1) ;
%     
%     SigmaL = sqrt(sum((LaplacianScale(:,i)-Fbar).^2)) / sqrt(NumPixels);
%     %  
%     LaplaceScaleInvariant(:,i) = ( LaplacianScale(:,i) - Fbar ) / SigmaL;
    
end


% Use this laplacian descriptor
LaplaceDOG2D = LaplacianScaleNormalized;


NumFeatures = 0;
FeaturePoint2D.Scale = zeros(10000,1);
FeaturePoint2D.Location = zeros(10000,1);

for i = 1 : NumPixels
    
    CurrentNeighbors = Neighbors2D{i,1};
    CurrentNeighbors = sub2ind([MaxImageSize(1),MaxImageSize(2)], CurrentNeighbors(:,1), CurrentNeighbors(:,2));
    
    for j = 2 : MaxLevel - 2
        
        CurrentValue = LaplaceDOG2D(i,j);
        
        SurroundingValues_under = LaplaceDOG2D([i;CurrentNeighbors],j-1);
        SurroundingValues_same = LaplaceDOG2D(CurrentNeighbors,j);
        SurroundingValues_above = LaplaceDOG2D([i;CurrentNeighbors],j+1);
        
        
        if all(CurrentValue > SurroundingValues_same)
            
            if all(CurrentValue > SurroundingValues_under)
                if all(CurrentValue > SurroundingValues_above)
                    NumFeatures = NumFeatures + 1;
                    % Maximum
                    FeaturePoint2D.Scale(NumFeatures, 1) = ScaleParamterImage(j);
                    FeaturePoint2D.Location(NumFeatures, 1) = i;
                end
            end
            
        elseif all(CurrentValue < SurroundingValues_same)
            
            if all(CurrentValue < SurroundingValues_under)
                if all(CurrentValue < SurroundingValues_above)
                    % Minimum
                    NumFeatures = NumFeatures + 1;
                    FeaturePoint2D.Scale(NumFeatures, 1) = ScaleParamterImage(j);
                    FeaturePoint2D.Location(NumFeatures, 1) = i;
                end
            end
        end
        
        
        % Check Both if it is maximum or minimum
        %         if all(CurrentValue > SurroundingValues)
        %             NumFeatures = NumFeatures + 1;
        %             % Maximum
        %             FeaturePoint2D.Scale(NumFeatures, 1) = ScaleParamterImage(j);
        %             FeaturePoint2D.Location(NumFeatures, 1) = i;
        %         end
        %         if all(CurrentValue < SurroundingValues)
        %             % Minimum
        %              NumFeatures = NumFeatures + 1;
        %             FeaturePoint2D.Scale(NumFeatures, 1) = ScaleParamterImage(j);
        %             FeaturePoint2D.Location(NumFeatures, 1) = i;
        %         end
        
    end
    
end

FeaturePoint2D.Scale(FeaturePoint2D.Location == 0) = [];
FeaturePoint2D.Location(FeaturePoint2D.Location == 0) = [];


% Comet67P Only
ZeroLogic = Image(FeaturePoint2D.Location) < 0.04;
FeaturePoint2D.Location(ZeroLogic,:) = [];
FeaturePoint2D.Scale(ZeroLogic,:) = [];

% save FeaturePointComet67P2DDiffusion2DConnectivitySigma FeaturePoint2D
% save Comet67P2DDiffusionSigma075 Signal2D

% save FeaturePointComet67P2D FeaturePoint2D
% save FeaturePointFlowers2D FeaturePoint2D
% load('FeaturePointFlowers2D.mat')
% load('FeaturePointComet67P2D.mat')

[SortValue2D, SortOrder2D] = sort(FeaturePoint2D.Scale,'descend');


CircleAngle = linspace(0, 2*pi, 360);
CircleX = cos(CircleAngle);
CircleY = sin(CircleAngle);




PlotFeatures = 1:500;%5:max(SortOrder2D);


figure
imshow(Image)
hold on



for i = 1 : length(PlotFeatures)
    
    [a,b] = ind2sub([MaxImageSize(1),MaxImageSize(2)], FeaturePoint2D.Location(SortOrder2D(i)));
    
    CircleAtPoint = bsxfun(@plus, SortValue2D(PlotFeatures(i))*[CircleX; CircleY], [a;b]);
    plot(CircleAtPoint(2,:), CircleAtPoint(1,:),'r')
% % % %     plot(CircleAtPoint(1,:), CircleAtPoint(2,:),'r')

    
    %     drawnow

end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Explicit 3D Surface, Laplace-Beltrami, diffusion for impulse 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NumPointSurf = length(1:eSS:NumVertex+(1/eSS -1)*eSS);
MiddleSurfPoint = round(NumPointSurf/2);

[xSurf3D, ySurf3D, zSurf3D] = meshgrid(1:eSS:MaxSurfSize(2)+(1/eSS -1)*eSS,1:eSS:MaxSurfSize(1)+(1/eSS -1)*eSS,0);


PointCloud.Location = [xSurf3D(:), ySurf3D(:), zSurf3D(:)];
% PointCloud.Location = [ySurf3D(:), xSurf3D(:), zSurf3D(:)];
PointCloud.LocationCount = length(PointCloud.Location);
% PointCloud.Face = delaunay(xSurf3D(:), ySurf3D(:));
PointCloud.Face = delaunay(ySurf3D(:), xSurf3D(:));
PointCloud.FaceCount = length(PointCloud.Face);
PointCloud.FaceArea = findFaceArea(PointCloud.Location, PointCloud.Face);
PointCloud.Signal = zeros(PointCloud.LocationCount,1);
PointCloud.Signal = reshape(Image, [],1);


Neighbors3D = findAdjacentNeighbors(PointCloud.Face, PointCloud.Location);

save_off(PointCloud.Location, PointCloud.Face, 'PlaneImage.off');

[LapMatMeshWeights, Area, hEdge2] = symmshlp_matrix('PlaneImage.off');

hEdge = (hEdge2/2);
% tauExplicit = 0.5*hEdge^2;
tauExplicit = hEdge/4;

MaxTauExplicit  = 100;
NumStepsExplcit = round(MaxTauExplicit / tauExplicit);

A1 = sparse(1:length(Area),1:length(Area), 1./Area);

LBM = A1 * LapMatMeshWeights;

ItL = speye(length(Area),length(Area)) - tauExplicit * LBM;

SignalExplicit = zeros(PointCloud.LocationCount, NumStepsExplcit);
SignalExplicit(:,1) = PointCloud.Signal;
SignalExplicitOld = PointCloud.Signal;

WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsExplcit-1));

for i = 1 : NumStepsExplcit - 1
    
    [SignalExplicitNew, flag] = bicg(ItL, SignalExplicitOld, 1e-10, 30);
    
    if flag
        flag
    end
    
    if any(SignalExplicitNew < 0) || any(SignalExplicitNew > 1)
        SignalExplicitNew(SignalExplicitNew < 0) = 0;
        SignalExplicitNew(SignalExplicitNew > 1) = 1;
        warning('Values out of bounds')
        i
    end
    
    SignalExplicit(:,i+1) = SignalExplicitNew;
    
    SignalExplicitOld = SignalExplicitNew;
    
    waitbar(i/NumStepsExplcit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsExplcit-1));
    
end

waitbar(i/NumStepsExplcit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)


% save SignalExplicitFlowers SignalExplicit


figure
for i = 1 : NumStepsExplcit
    a = reshape(SignalExplicit(:,i),MaxSurfSize(1),MaxSurfSize(2));

    imshow(a)
    
end



% Find the scale parameters
maxsample = 2; 
ws = 0 : 0.001 : maxsample;
NumSample = length(ws);
ScaleParameterSpatialExplicit = zeros(NumStepsExplcit,2);
cut = sqrt(log(2));
db3 = 1/sqrt(2);
ws2 = (ws.^2)';
H =  ones(NumSample,1) ;
h = 1 ./ ( ones(NumSample,1) + tauExplicit * ws2 );

for i = 2 : NumStepsExplcit
    
    % Transfer function
    % 1st order
    H = H .* h;
    % Find the frequency at the cutoff values
    CutoffFrequency = interp1(H, ws, db3);
    % Change cutoff frequency to scale parameter
    ScaleParameterFrequency = CutoffFrequency / cut;  
    ScaleParameterSpatialExplicit(i,1) = 1 / ScaleParameterFrequency;
    ScaleParameterSpatialExplicit(i,2) = sqrt(2*i*tauExplicit);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimated Laplacians (Difference of Gaussian DoG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defined by Fadaifard
% Scale Normalized Laplacian, and scale Invariant Laplacian

LaplacianScale = zeros(PointCloud.LocationCount, NumStepsExplcit-1);
LaplacianScaleNormalized = zeros(PointCloud.LocationCount, NumStepsExplcit-1);
LaplaceScaleInvariant = zeros(PointCloud.LocationCount, NumStepsExplcit-1);


for i = 1 : NumStepsExplcit-1
    
    % Estimated Laplacian
    LaplacianScale(:,i) = ( SignalExplicit(:,i+1) - SignalExplicit(:,i) ) / (ScaleParameterSpatialExplicit(i+1,1)^2-ScaleParameterSpatialExplicit(i,1)^2);
    
    % Scale Normalized Laplacian
    LaplacianScaleNormalized(:,i) = 2 * ( SignalExplicit(:,i+1) - SignalExplicit(:,i) ) * (ScaleParameterSpatialExplicit(i,1)^2 / (ScaleParameterSpatialExplicit(i+1,1)^2 - ScaleParameterSpatialExplicit(i,1)^2) );
    
    % Scale Invariant Laplace
    Fbar = (sum(SignalExplicit(:,i)) / PointCloud.LocationCount) * ones(PointCloud.LocationCount, 1) ;
    
    SigmaL = sqrt(sum((LaplacianScale(:,i)-Fbar).^2)) / sqrt(PointCloud.LocationCount);
    %  
    LaplaceScaleInvariant(:,i) = ( LaplacianScale(:,i) - Fbar ) / SigmaL;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extrema Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this laplacian descriptor
LaplaceDOG = LaplacianScaleNormalized;

FeaturePoint.Scale = zeros(11000,2);
FeaturePoint.Location = zeros(11000,1);

NumFeatures = 0;

for i = 1 : PointCloud.LocationCount
    
    CurrentNeighbors = Neighbors3D{i,1};
    
    for j = 2 : NumStepsExplcit - 2
        
        CurrentValue = LaplaceDOG(i,j);
        
        SurroundingValues = [LaplaceDOG([i,CurrentNeighbors],j-1);
            LaplaceDOG(CurrentNeighbors,j);
            LaplaceDOG([i,CurrentNeighbors],j+1)];
        
        % Check Both if it is maximum or minimum
        if all(CurrentValue > SurroundingValues)
            NumFeatures = NumFeatures + 1;
            % Maximum
            FeaturePoint.Scale(NumFeatures, 1:2) = ScaleParameterSpatialExplicit(j,:);
            FeaturePoint.Location(NumFeatures, 1) = i;
        end
        if all(CurrentValue < SurroundingValues)
            % Minimum
             NumFeatures = NumFeatures + 1;
            FeaturePoint.Scale(NumFeatures, 1:2) = ScaleParameterSpatialExplicit(j,:);
            FeaturePoint.Location(NumFeatures, 1) = i;
        end
        
    end
    
end



FeaturePoint.Scale(FeaturePoint.Location == 0,:) = [];
FeaturePoint.Location(FeaturePoint.Location == 0) = [];

% save FeaturePointFlowersExplicit3D FeaturePoint



[SortValue3D, SortOrder3D] = sort(FeaturePoint.Scale(:,2),'descend');


CircleAngle = linspace(0, 2*pi, 360);
CircleX = cos(CircleAngle);
CircleY = sin(CircleAngle);



PlotFeatures = 1:200;%:max(SortOrder3D);


figure
imshow(reshape(PointCloud.Signal,MaxSurfSize,MaxSurfSize))
hold on



for i = 1 : length(PlotFeatures)
% % % % %     CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY], [PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),1);PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),2)]);
    CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY], [PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),2);PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),1)]);

    plot(CircleAtPoint(1,:), CircleAtPoint(2,:),'r')
    

%     pause
end


for i = [22, 44, 86]
   figure
   imshow(reshape(SignalExplicit(:,i),MaxSurfSize,MaxSurfSize))
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Implicit 3D Surface, L, E, M, diffusion for impulse 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create 3D implcit Surface
bandwidth = 1.00001*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));

MinPoint = round(min(PointCloud.Location) - bandwidth - spacing, 1);
MaxPoint = round(max(PointCloud.Location) + bandwidth + spacing, 1);

x1d = (MinPoint(1):spacing:MaxPoint(1))';
y1d = (MinPoint(2):spacing:MaxPoint(2))';
z1d = (MinPoint(3):spacing:MaxPoint(3))';

[IJK,DIST,CP,XYZ,CPFACE] = tri2cp(PointCloud.Face, PointCloud.Location, spacing, MinPoint, porder, 1);

% MinXYZ = min(XYZ);
% MaxXYZ = max(XYZ);
% MinIJK = min(IJK);
% MaxIJK = max(IJK);



% XYZ = MinPoint - spacing + IJK * spacing

BandSearchSize = [length(x1d), length(y1d), length(z1d)];

Band = sub2ind(BandSearchSize, IJK(:,1), IJK(:,2), IJK(:,3));


% Create L, E, M
L = laplacian_3d_matrix(y1d,x1d,z1d, Lorder, Band);

% Eplot = interpLagrange3D(BandSearchSize, MinPoint, PointCloud.Location(:,[1,2,3]), porder, Band, spacing);

Eplot = interp3_matrix(y1d, x1d, z1d, PointCloud.Location(:,2), PointCloud.Location(:,1), PointCloud.Location(:,3), porder, Band);


% Ecp = interpLagrange3D(BandSearchSize, MinPoint, CP(:,[2,1,3]), porder, Band, spacing);
Ecp = interp3_matrix(y1d, x1d, z1d, CP(:,2), CP(:,1), CP(:,3), porder, Band);


% save Ecp67PFlip Ecp
% save Eplot67PFlip Eplot
% save L67PFlip L
% save IJK67PFlip IJK
% save CP67PFlip CP
% save CPFACE67PFlip CPFACE

% load('Eplot67P.mat')
% load('CP67P.mat')
% load('IJK67P.mat')
% load('CPFACE67P.mat')
% load('L67P.mat')
% load('Ecp67P.mat')

load('Eplot67PFlip.mat')
load('CP67PFlip.mat')
load('IJK67PFlip.mat')
load('CPFACE67PFlip.mat')
load('L67PFlip.mat')
load('Ecp67PFlip.mat')

M = lapsharp(L, Ecp);


% Extrapolate data to embedding

FaceInterpolateWeights = interpBarycenterTriangle(PointCloud, CP, CPFACE);

Signal = FaceInterpolateWeights * PointCloud.Signal;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precompute the scale parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the scale parameters
% maxsample = 10; 
% ws = 0 : 0.01 : maxsample;
% NumSample = length(ws);
ScaleParameterSpatialImplicit = zeros(NumStepsImplicit,1);
% cut = sqrt(log(2));
% db3 = 1/sqrt(2);
% ws2 = (ws.^2)';
% H = ones(NumSample,1) ;
% h = 1 ./ ( ones(NumSample,1) + tauImplicit/spacing^2 * ws2 );

for i = 1 : NumStepsImplicit - 1
    
    % Transfer function
    % 1st order
%     H = H .* h;
    % Find the frequency at the cutoff values
%     [uH, indH] = unique(H);
%     CutoffFrequencyImplicit = interp1(uH, ws(indH), db3);
    % Change cutoff frequency to scale parameter
%     ScaleParameterFrequencyImplicit = CutoffFrequencyImplicit / cut;   
%     ScaleParameterSpatialImplicit(i+1,1) =  spacing / ScaleParameterFrequencyImplicit;
    ScaleParameterSpatialImplicit(i+1,1) = sqrt(2*i*tauImplicit);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SignalImplicit = zeros(PointCloud.LocationCount, NumStepsImplicit);
SignalImplicit(:,1) = PointCloud.Signal;


ItM = speye(size(M)) - tauImplicit * M;

WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsImplicit-1));

for i = 1 : NumStepsImplicit - 1
    

    [SignalNew, flag, relres] = bicg(ItM, Signal, 1e-10, 30);

%     [SignalNew, flag, relres] = gmres(ItM, Signal, [], 1e-10, 30);

    if flag
        flag
        relres
    end
    
    if any(SignalNew < 0) || any(SignalNew > 1) 
        SignalNew(SignalNew < 0) = 0;
        SignalNew(SignalNew > 1) = 1;
        warning('Values out of bounds')
        i
    end
    
    
    Signal = SignalNew;
    
    SignalImplicit(:,i+1) = Eplot * SignalNew;
    waitbar(i/NumStepsImplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsImplicit-1));
    
end

waitbar(i/NumStepsImplicit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)


if any(SignalImplicit(:) < 0) || any(SignalImplicit(:) > 1) 
    SignalImplicit(SignalImplicit < 0) = 0;
    SignalImplicit(SignalImplicit > 1) = 1;
    warning('Values out of bounds')
end

save Comet67P3DDiffusionSigma075Flipped SignalImplicit

figure
for i = 1 : NumStepsImplicit
    a = reshape(SignalImplicit(:,i),MaxSurfSize(1),MaxSurfSize(2));

    imshow(a)
%     pause
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimated Laplacians (Difference of Gaussian DoG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defined by Fadaifard
% Scale Normalized Laplacian, and scale Invariant Laplacian

LaplacianScale = zeros(PointCloud.LocationCount, NumStepsImplicit-1);
LaplacianScaleNormalized = zeros(PointCloud.LocationCount, NumStepsImplicit-1);
LaplaceScaleInvariant = zeros(PointCloud.LocationCount, NumStepsImplicit-1);


for i = 1 : NumStepsImplicit-1
    
    % Estimated Laplacian
    LaplacianScale(:,i) = ( SignalImplicit(:,i+1) - SignalImplicit(:,i) ) / (ScaleParameterSpatialImplicit(i+1,1)^2-ScaleParameterSpatialImplicit(i,1)^2);
    
    % Scale Normalized Laplacian
    LaplacianScaleNormalized(:,i) = 2 * ( SignalImplicit(:,i+1) - SignalImplicit(:,i) ) * (ScaleParameterSpatialImplicit(i,1)^2 / (ScaleParameterSpatialImplicit(i+1,1)^2 - ScaleParameterSpatialImplicit(i,1)^2) );
    
    % Scale Invariant Laplace
    Fbar = (sum(SignalImplicit(:,i)) / PointCloud.LocationCount) * ones(PointCloud.LocationCount, 1) ;
    
    SigmaL = sqrt(sum((LaplacianScale(:,i)-Fbar).^2)) / sqrt(PointCloud.LocationCount);
    %  
    LaplaceScaleInvariant(:,i) = ( LaplacianScale(:,i) - Fbar ) / SigmaL;
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extrema Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this laplacian descriptor
LaplaceDOG = LaplacianScaleNormalized;

FeaturePoint.Scale = zeros(10000,2);
FeaturePoint.Location = zeros(10000,1);

NumFeatures = 0;

for i = 1 : PointCloud.LocationCount
    
    CurrentNeighbors = Neighbors3D{i,1};
    
    for j = 2 : NumStepsImplicit - 2
        
        CurrentValue = LaplaceDOG(i,j);
        
        SurroundingValues_under = LaplaceDOG([i,CurrentNeighbors],j-1);
        SurroundingValues_same = LaplaceDOG(CurrentNeighbors,j);
        SurroundingValues_above = LaplaceDOG([i,CurrentNeighbors],j+1);
        
        
        if all(CurrentValue > SurroundingValues_same)
            
            if all(CurrentValue > SurroundingValues_under)
                if all(CurrentValue > SurroundingValues_above)
                    NumFeatures = NumFeatures + 1;
                    % Maximum
                    FeaturePoint.Scale(NumFeatures, 1:2) = ScaleParameterSpatialImplicit(j,:);
                    FeaturePoint.Location(NumFeatures, 1) = i;
                end
            end
            
        elseif all(CurrentValue < SurroundingValues_same)
            
            if all(CurrentValue < SurroundingValues_under)
                if all(CurrentValue < SurroundingValues_above)
                    % Minimum
                    NumFeatures = NumFeatures + 1;
                    FeaturePoint.Scale(NumFeatures, 1:2) = ScaleParameterSpatialImplicit(j,:);
                    FeaturePoint.Location(NumFeatures, 1) = i;
                end
            end
        end
        
        
        % Check Both if it is maximum or minimum
%         if all(CurrentValue > SurroundingValues)
%             NumFeatures = NumFeatures + 1;
%             % Maximum
%             FeaturePoint.Scale(NumFeatures, 1:2) = ScaleParameterSpatialImplicit(j,:);
%             FeaturePoint.Location(NumFeatures, 1) = i;
%         end
%         if all(CurrentValue < SurroundingValues)
%             % Minimum
%              NumFeatures = NumFeatures + 1;
%             FeaturePoint.Scale(NumFeatures, 1:2) = ScaleParameterSpatialImplicit(j,:);
%             FeaturePoint.Location(NumFeatures, 1) = i;
%         end
        
    end
    
end

FeaturePoint.Scale(FeaturePoint.Location == 0,:) = [];
FeaturePoint.Location(FeaturePoint.Location == 0) = [];

% For Comet67P
% LowLightLogic = Image(FeaturePoint.Location)< 0.04;
% FeaturePoint.Location(LowLightLogic,:) = [];
% FeaturePoint.Scale(LowLightLogic,:) = [];

% save FeaturePointComet67PImplicit3D FeaturePoint

% load('FeaturePointCameraMan.mat')
% load('SignalImplicitCameraMan.mat')

% load('SignalImplicitFlowers.mat')
% load('FeaturePointFlowers3D.mat')

% load('SignalImplicitCome67P.mat')
% load('FeaturePointComent67PImplicit3D.mat')

load('FeaturePointFlowers_FullRes.mat')

[SortValue3D, SortOrder3D] = sort(FeaturePoint.Scale(:,2),'descend');


CircleAngle = linspace(0, 2*pi, 360);
CircleX = cos(CircleAngle);
CircleY = sin(CircleAngle);




PlotFeatures = 1:500;%:max(SortOrder3D);


figure
imshow(reshape(PointCloud.Signal,MaxSurfSize(1),MaxSurfSize(2)))
hold on



for i = 1 : length(PlotFeatures)
    CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY], [PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),1);PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),2)]);
% % % % % % %     CircleAtPoint = bsxfun(@plus, SortValue3D(PlotFeatures(i))*[CircleX; CircleY], [PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),2);PointCloud.Location(FeaturePoint.Location(SortOrder3D(PlotFeatures(i))),1)]);

    plot(CircleAtPoint(1,:), CircleAtPoint(2,:),'r')
    

%     pause
end



for i = [10, 50, 100, 200, NumStepsImplicit]
   figure
   imshow(reshape(SignalImplicit(:,i),MaxSurfSize(1),MaxSurfSize(2)))
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matching Key Points from 2D / 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






FeaturePoint
FeaturePoint2D


for i = 1 : length(FeaturePoint2D.Scale)
    
    [a,b] = ind2sub([MaxImageSize, MaxImageSize], FeaturePoint2D.Location(i));

    for j = 1 : length(FeaturePoint.Scale);
        
        LocationDistance = bsxfun(@minus, FeaturePoint.Location, [a,b]);
        
        
        
    end
    
    
end






























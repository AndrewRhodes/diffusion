


close all
clear
clc

global ProjectRoot; % Additional Paths


MaxImageSize = 100;
MiddleImage = MaxImageSize / 2;
NumSteps = 300;
sigma2D = 0.75;

ScaleParameter2D = findScaleParameter(sigma2D, 1, NumSteps, 'Gaussian', 'Natural');

gausswindow = 5;

Gauss = fspecial('gaussian', [gausswindow, gausswindow], sigma2D);



[xImage2D, yImage2D] = ndgrid(1:MaxImageSize,1:MaxImageSize);

xyImage2D = [xImage2D(:), yImage2D(:)];

Image2Dcenter = sub2ind([MaxImageSize, MaxImageSize], MiddleImage, MiddleImage);

RadialDist2D = sqrt(sum(bsxfun(@minus, xyImage2D, xyImage2D(Image2Dcenter,:)).^2,2));




Signal2D = zeros(MaxImageSize, MaxImageSize, NumSteps);
Signal2D(MiddleImage, MiddleImage, 1) = 1;


figure
for i = 1 : NumSteps - 1
   Signal2D(:,:,i+1) = imfilter(Signal2D(:,:,i), Gauss, 'replicate', 'same', 'conv');
%    Signal2D(:,:,i+1) = Signal2D(:,:,i+1) ./ max(max(Signal2D(:,:,i+1)));
%    imshow(Signal2D(:,:,i+1),[])
end



figure
xgauss = 0:0.01:20;
for i = 1 : NumSteps -1
    clf
    
    Gaus = exp( -(xgauss.^2) / (2*ScaleParameter2D(i)^2) );
    Gaus = Gaus * max(max(Signal2D(:,:,i)));
    
    plot(xgauss, Gaus, 'b:') 
    hold on
    plot(RadialDist2D, reshape(Signal2D(:,:,i),[],1),'r.')
    xlim([0 10])
    
    ylabel('Intensity')
    xlabel('Radial Distance')

    pause(0.1)
end













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumSteps = 100;
alpha = 1;
BDF = 1;

options.rho = 6;
options.dtype = 'euclidean'; % 'euclidean', 'geodesic' %
options.htype = 'psp'; % 'psp', 'ddr'

setTau = @(e_bar) e_bar^(2/5);
setHs = @(e_bar) 2*e_bar^(1/5);

Model = strcat('Plane_',num2str(MaxImageSize),'x',num2str(MaxImageSize));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eSS = 0.5;
MaxSurfSize = 100;
NumPointSurf = length(1:eSS:MaxSurfSize+(1/eSS -1)*eSS);
MiddleSurfPoint = round(NumPointSurf/2);
PointCloudCenter = sub2ind([NumPointSurf, NumPointSurf], MiddleSurfPoint, MiddleSurfPoint);


FileNameOff = strcat(ProjectRoot,'/models/',Model,'_e',num2str(eSS),'.off');

if ~exist(FileNameOff, 'file')
    [xSurf3D, ySurf3D, zSurf3D] = meshgrid(1:eSS:MaxSurfSize+(1/eSS -1)*eSS,1:eSS:MaxSurfSize+(1/eSS -1)*eSS,0);
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
PointCloud.Signal = zeros(PointCloud.LocationCount,1);
PointCloud.Signal(PointCloudCenter) = 1;
PointCloud = findMeshResolution(PointCloud, 'Model');


% FileNameNeighbors = strcat(ProjectRoot, '/models/Neighbors_',Model,'.mat');

% if ~exist( FileNameNeighbors, 'file')
%     [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);
%     save(FileNameNeighbors ,'Neighbors', '-v7.3')
% else
%     load(FileNameNeighbors, 'Neighbors');
% end

% PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);

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
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[LapMatMeshWeights, Area, hEdge2] = symmshlp_matrix( FileNameOff, options);

Alength = length(Area);

A1 = sparse(1:Alength, 1:Alength, 1./Area);

LBM1 = A1 * LapMatMeshWeights;

LBMdiag = abs(diag(1./diag(LBM1)));

LBM2 = LBMdiag * LBM1;

ItL1{1,1} = speye(Alength, Alength) - alpha * tau * LBM1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



ScaleParameter3D = findScaleParameter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

Signal3D = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);



figure
for i = 1 : NumSteps - 1
%    Signal2D(:,:,i+1) = imfilter(Signal2D(:,:,i), Gauss, 'replicate', 'same', 'conv');
%    Signal2D(:,:,i+1) = Signal2D(:,:,i+1) ./ max(max(Signal2D(:,:,i+1)));
%    imshow(Signal2D(:,:,i+1),[])
   imshow(reshape(Signal3D(:,i+1), NumPointSurf, NumPointSurf), [])
end



RadialDist3D = sqrt(sum(bsxfun(@minus, PointCloud.Location, PointCloud.Location(PointCloudCenter,:)).^2,2));
figure
xgauss = 0:0.01:20;
for i = 1 : NumSteps -1
    clf
    
    Gaus = exp( -(xgauss.^2) / (2*ScaleParameter3D(i)^2) );
    Gaus = Gaus * max(max(Signal3D(:,i)));
    
    plot(xgauss, Gaus, 'b:') 
    hold on
    plot(RadialDist3D, Signal3D(:,i),'r.')
    xlim([0 20])
    
    ylabel('Intensity')
    xlabel('Radial Distance')
    title(sprintf('Iter: %d, Scale: %0.3f',i, ScaleParameter3D(i)))

    pause(0.5)
end


















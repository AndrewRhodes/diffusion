% Andrew Rhodes
% ASEL
% 2018


close all
clear 
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional Paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd ~/Desktop/Ashish/'CS3 Code'/
addpath(genpath('~/Documents/Software/MeshLP/'))
addpath(genpath('~/Documents/Software/cp_matrices/'))
addpath('~/AFOSR/Ashish/CS3 Code/')
addpath(genpath('~/GitProjects/pose/MATLAB_PointCloudDescriptors/OURCVFH/'))
addpath('src/')
addpath('data/')
addpath('models/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NumSteps = 500;

MaxImageSize = [30 30];
MiddleImage = MaxImageSize / 2;


eSS = 0.5; % Explicit Surface Spacing
MaxSurfSize = [30 30];




NumPointSurf = length(1:eSS:MaxSurfSize(1)+(1/eSS -1)*eSS);
MiddleSurfPoint = round(NumPointSurf/2);

[xSurf3D, ySurf3D, zSurf3D] = ndgrid(1:eSS:MaxSurfSize(1)+(1/eSS -1)*eSS,1:eSS:MaxSurfSize(2)+(1/eSS -1)*eSS,0);

PointCloudCenter = sub2ind([NumPointSurf, NumPointSurf], MiddleSurfPoint, MiddleSurfPoint);

PointCloud.Location = [xSurf3D(:), ySurf3D(:), zSurf3D(:)];
PointCloud.LocationCount = length(PointCloud.Location);
PointCloud.Face = delaunay(xSurf3D(:), ySurf3D(:));
PointCloud.FaceCount = length(PointCloud.Face);
PointCloud.FaceArea = findFaceArea(PointCloud.Location, PointCloud.Face);
PointCloud.Signal = zeros(PointCloud.LocationCount,1);
PointCloud.Signal(PointCloudCenter,1) = 1;
PointCloud.VertexArea = accumarray( reshape(PointCloud.Face, 3*PointCloud.FaceCount ,1), repmat(PointCloud.FaceArea,3,1), [PointCloud.LocationCount, 1] )/3;
PointCloud = findMeshResolution(PointCloud, 'Model');


AverageEdgeLength = findAverageEdgeLengths(PointCloud.Location, PointCloud.Face, 1);

tic
G = makeMeshGaussian(PointCloud, AverageEdgeLength, 6*sqrt(2)*AverageEdgeLength, 1);
toc


 
SignalExplicit = zeros(PointCloud.LocationCount, NumSteps);
SignalExplicit(:,1) = PointCloud.Signal;



WaitBar = waitbar(0, sprintf('Exp Gaus Diff Step %i of %i', 0, NumSteps-1));
figure
for i = 1 : NumSteps - 1
    
    SignalExplicit(:,i+1) = G * SignalExplicit(:,i);
    
    if flag
        flag
    end
    
    
     imshow(reshape(SignalExplicit(:,i), NumPointSurf, NumPointSurf),[])
     drawnow
     
    waitbar(i/NumSteps, WaitBar, sprintf('Exp Gauss Diff Step %i of %i', i, NumSteps-1));
    
end
waitbar(i/PointCloud.LocationCount, WaitBar, sprintf('Exp Gauss Diff Step Complete'));
close(WaitBar)












RadialDist3D = sqrt(sum(bsxfun(@minus, PointCloud.Location, PointCloud.Location(PointCloudCenter,:)).^2,2));

xgauss = 0:0.01:20;
figure
for i = 1 : NumSteps
    
    
    clf
    plot(RadialDist3D, SignalExplicit(:,i),'r.')
    hold on
    axis([0 20 0 2*max(SignalExplicit(:,i))])
    
    t = (i-1) * AverageEdgeLength^2 * median(PointCloud.VertexArea)^(1/3)
%     t = (i-1) * ((AverageEdgeLength)/median(PointCloud.VertexArea)^(2/3))^2  ; %((AverageEdgeLength^(2/3))*median(PointCloud.VertexArea)^(1/3))

    sqrt(t)
%     t = 4
    
    Gauss = exp(-(xgauss.^2) / (2*t));
    Gauss = Gauss * max(SignalExplicit(:,i));
    plot(xgauss, Gauss,'b-')
    
    pause
end










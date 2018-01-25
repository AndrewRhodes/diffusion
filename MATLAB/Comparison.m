% Andrew Rhodes
% ASEL
% 2017

% Perform scale space on an 2D image of an impulse and of a scene.
% Sample the image onto a 3D explicit surface.
% Perform explicit surface diffusion and compare.
% Perform implicit surface diffusion and copare.


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

% 2D image diffusion
MaxLevel = 70;
gausswindow = 2;
tau2D = 0.75;

MaxImageSize = 30;
MiddleImage = MaxImageSize / 2;


% 3D explcit surface diffusion
eSS = 1; % Explicit Surface Spacing
MaxSurfSize = 30;


% 3D implicit surface
spacing = 0.25;
sigma = 0.25;
porder = 3;
Lorder = 2;
dim = 3;
tauImplicit = spacing / 4;
tauImplicit = tau2D^2/2;
MaxTauImplicit = 70;
NumStepsImplicit = round(MaxTauImplicit / tauImplicit);
% NumStepsImplicit = MaxTauImplicit / sigma^2


% If using Gaussian implicit surface diffusion, p=3 causes negative values
% porder = 4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Impulse Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Signal2D = zeros(MaxImageSize, MaxImageSize, MaxLevel);
Signal2D(MiddleImage, MiddleImage, 1) = 1;

% Gauss = fspecial('gaussian', [gausswindow, gausswindow], tau2D);


[x,y] = meshgrid(-gausswindow:gausswindow,-gausswindow:gausswindow);
Gauss =  (1/(2*pi*tau2D^2)) .* exp( -(x.^2 + y.^2) ./ (2*tau2D^2) );
Gauss = Gauss./sum(Gauss(:))

for i = 1 : MaxLevel - 1
    
%     imshow(Signal2D(:,:,i),[])
%     pause
    
%     Signal2D(:,:,i+1) = imgaussfilt(Signal2D(:,:,1), sqrt(i)*tau2D);
%     Signal2D(:,:,i+1) = imgaussfilt(Signal2D(:,:,i), tau2D);
    Signal2D(:,:,i+1) = imfilter(Signal2D(:,:,i), Gauss,'replicate');

    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 2D Diffusion at Scales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xImage2D, yImage2D] = ndgrid(1:MaxImageSize,1:MaxImageSize);

xyImage2D = [xImage2D(:), yImage2D(:)];

Image2Dcenter = sub2ind([MaxImageSize, MaxImageSize], MiddleImage, MiddleImage);

RadialDist2D = sqrt(sum(bsxfun(@minus, xyImage2D, xyImage2D(Image2Dcenter,:)).^2,2));




%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures for paper %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
xgauss = 0:0.01:20;
for i = 1 : MaxLevel -1
    clf
    i
    t2 = (i-1) * tau2D^2
    
    Gaus = exp( -(xgauss.^2) / (2*t2) );
    Gaus = Gaus  * max(max(Signal2D(:,:,i)));
    
    plot(xgauss, Gaus, 'b:') 


    hold on
   
    plot(RadialDist2D, reshape(Signal2D(:,:,i),[],1),'r.')
    hold on
    xlim([0 10])
    
    ylabel('Intensity')
    xlabel('Radial Distance')

    pause
end

% % tau2D = 0.5; [6,21,61]
% % tau2d = 0.75; [5]

xgauss = 0:0.01:20;
for i = [5,15,35]
%     clf
    figure

    i;
    t2 = (i-1) * tau2D^2
    
    Gaus = exp( -(xgauss.^2) / (2*t2) );
    Gaus = Gaus  * max(max(Signal2D(:,:,i)));
    
    plot(xgauss, Gaus, 'b:','linewidth',3);

    hold on
    
    plot(RadialDist2D, reshape(Signal2D(:,:,i),[],1), 'r.','markersize',30)

    hold on
    xlim([0 10])
    
    ylabel('Intensity')
    xlabel('Radial Distance')
    ax = gca;
    ax.XAxis.FontSize = 55;
    ax.YAxis.FontSize = 55;
%     pause
end


figure
imshow(Signal2D(:,:,5),[])

figure
imshow(Signal2D(:,:,15),[])

figure 
imshow(Signal2D(:,:,35),[])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Explicit 3D Surface, Laplace-Beltrami, diffusion for impulse 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NumPointSurf = length(1:eSS:MaxSurfSize+(1/eSS -1)*eSS);
MiddleSurfPoint = round(NumPointSurf/2);

[xSurf3D, ySurf3D, zSurf3D] = ndgrid(1:eSS:MaxSurfSize+(1/eSS -1)*eSS,1:eSS:MaxSurfSize+(1/eSS -1)*eSS,0);

PointCloudCenter = sub2ind([NumPointSurf, NumPointSurf], MiddleSurfPoint, MiddleSurfPoint);

PointCloud.Location = [xSurf3D(:), ySurf3D(:), zSurf3D(:)];
PointCloud.LocationCount = length(PointCloud.Location);
PointCloud.Face = delaunay(xSurf3D(:), ySurf3D(:));
PointCloud.FaceCount = length(PointCloud.Face);
PointCloud.FaceArea = findFaceArea(PointCloud.Location, PointCloud.Face);
PointCloud.Signal = zeros(PointCloud.LocationCount,1);
PointCloud.Signal(PointCloudCenter,1) = 1;
PointCloud = findMeshResolution(PointCloud, 'Model');

save_off(PointCloud.Location, PointCloud.Face, strcat('Plane',num2str(NumPointSurf),'.off'));

[LapMatMeshWeights, Area, hEdge2] = symmshlp_matrix(strcat('Plane',num2str(NumPointSurf),'.off'));

% [LapMatMeshWeights, Area, hEdge2] = symmshlp_matrix('Plane120.off');
hEdge = (hEdge2/2);
% tauExplicit = 0.5*hEdge^2;
tauExplicit = hEdge/4;

MaxTauExplicit  = 40;
NumStepsExplcit = round(MaxTauExplicit / tauExplicit);

A1 = sparse(1:length(Area),1:length(Area), 1./Area);

LBM = A1 * LapMatMeshWeights;

ItL = speye(length(Area),length(Area)) - tauExplicit * LBM;

SignalExplicit = zeros(PointCloud.LocationCount, NumStepsExplcit);
SignalExplicit(:,1) = PointCloud.Signal;


WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsExplcit-1));

for i = 1 : NumStepsExplcit - 1
    
    [SignalExplicit(:,i+1), flag] = bicg(ItL, SignalExplicit(:,i), 1e-10, 30);
    if flag
        flag
    end
    
    waitbar(i/NumStepsExplcit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsExplcit-1));
    
end

waitbar(i/NumStepsExplcit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)



for i = 1 : NumStepsExplcit 
    
    imshow(reshape(SignalExplicit(:,i), NumPointSurf, NumPointSurf),[])
 
    
end




% Find the scale parameters
maxsample = 2; 
ws = 0 : 0.001 : maxsample;
NumSample = length(ws);
ScaleParameterSpatial = zeros(NumStepsExplcit,1);
ScaleParameterSpatial2 = zeros(NumStepsExplcit,1);
cut = sqrt(log(2));
db3 = 1/sqrt(2);
ws2 = (ws.^2)';
H =  ones(NumSample,1) ;
H2 =  ones(NumSample,1) ;
h = 1 ./ ( ones(NumSample,1) + tauExplicit/hEdge^2 * ws2 );
h2 = 1 ./ ( ones(NumSample,1) + tauExplicit * ws2 );

for i = 2 : NumStepsExplcit
    
    % Transfer function
    % 1st order
    H = H .* h;
    H2 = H2 .* h2;
    % Find the frequency at the cutoff values
    CutoffFrequency = interp1(H, ws, db3);
    % Change cutoff frequency to scale parameter
    ScaleParameterFrequency = CutoffFrequency / cut;   
    ScaleParameterSpatial(i,1) = hEdge / ScaleParameterFrequency;


    CutoffFrequency = interp1(H2, ws, db3);
    % Change cutoff frequency to scale parameter
    ScaleParameterFrequency = CutoffFrequency / cut;  
    ScaleParameterSpatial2(i,1) = 1 / ScaleParameterFrequency;

end




% tauExplicit = (hEdge/2)^2;
% ScaleParameterSpatial(i,1)
% sqrt(2*tauExplicit*i)


RadialDist3D = sqrt(sum(bsxfun(@minus, PointCloud.Location, PointCloud.Location(PointCloudCenter,:)).^2,2));

xgauss = 0:0.01:20;
for i = 1 : NumStepsExplcit
    clf
    plot(RadialDist3D, SignalExplicit(:,i),'r.')
    hold on
    axis([0 10 0 2*max(SignalExplicit(:,i))])
    Gauss = exp(-(xgauss.^2) / (2*ScaleParameterSpatial(i)^2));
    Gauss = Gauss * max(SignalExplicit(:,i));
    plot(xgauss, Gauss,'b-')
    
    Gauss2 = exp(-(xgauss.^2) / (2*ScaleParameterSpatial2(i)^2));
    Gauss2 = Gauss2 * max(SignalExplicit(:,i));
    plot(xgauss, Gauss2,'g:')
    
    i
    ScaleParameterSpatial(i)^2
    ScaleParameterSpatial2(i)^2
    t = (i-1)*2*tauExplicit
   
    Gauss2 = exp(-(xgauss.^2) / (2*t));
    Gauss2 = Gauss2 * max(SignalExplicit(:,i));
    plot(xgauss, Gauss2,'k:')
    
    pause
end

% scale^2 = [8 16 35]
% eSS = 1, t=h/4 [15 29 63]
% eSS = 0.5, t=h/4 [29 57 124]
% eSS = 0.25 t=h/4 [58 113 247]

% [7 23 67]

% [3 10 28]
%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures for paper %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RadialDist3D = sqrt(sum(bsxfun(@minus, PointCloud.Location, PointCloud.Location(PointCloudCenter,:)).^2,2));

xgauss = 0:0.01:20;
% eSS = 1, t=h/4 [15 29 63]
% eSS = 0.5, t=h/4 [29 57 124]
% eSS = 0.25 t=h/4 [58 113 247]
for i = [58 113 247]
    
    figure
    plot(RadialDist3D, SignalExplicit(:,i),'r.','markersize',30)
    hold on
%     axis([0 10 0 2*max(SignalExplicit(:,i))])
    Gauss = exp(-(xgauss.^2) / (2*ScaleParameterSpatial(i)^2));
    Gauss = Gauss * max(SignalExplicit(:,i));
    plot(xgauss, Gauss,'m--','linewidth',3)
    
    ScaleParameterSpatial(i)^2
    t = (i-1)*2*tauExplicit
   
    Gauss2 = exp(-(xgauss.^2) / (2*t));
    Gauss2 = Gauss2 * max(SignalExplicit(:,i));
    plot(xgauss, Gauss2,'b:','linewidth',3)
    xlim([0 10])
    ylabel('Intensity')
    xlabel('Radial Distance')
    ax = gca;
    ax.XAxis.FontSize = 55;
    ax.YAxis.FontSize = 55;
%     pause
end


% eSS = 1, t=h/4 [15 29 63]
% eSS = 0.5, t=h/4 [29 57 124]
% eSS = 0.25 t=h/4 [58 113 247]
for i = [58 113 247]
    
    figure
    imshow(reshape(SignalExplicit(:,i), NumPointSurf, NumPointSurf),[])

    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Implicit 3D Surface, L, E, M, diffusion for impulse 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create 3D implcit Surface
bandwidth = 1.00001*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));
% Implicit Clostest Point Method
MinPoint = round(min(PointCloud.Location) - bandwidth - spacing, 1);
MaxPoint = round(max(PointCloud.Location) + bandwidth + spacing, 1);

% Gaussian Method
% MinPoint = round(min(PointCloud.Location) - bandwidth - 2*spacing, 1);
% MaxPoint = round(max(PointCloud.Location) + bandwidth + 2*spacing, 1);

x1d = (MinPoint(1):spacing:MaxPoint(1))';
y1d = (MinPoint(2):spacing:MaxPoint(2))';
z1d = (MinPoint(3):spacing:MaxPoint(3))';

[IJK,DIST,CP,XYZ,CPFACE] = tri2cp(PointCloud.Face, PointCloud.Location, spacing, MinPoint, porder, Lorder/2);

% MinXYZ = min(XYZ);
% MaxXYZ = max(XYZ);
% MinIJK = min(IJK);
% MaxIJK = max(IJK);



% XYZ = MinPoint - spacing + IJK * spacing

BandSearchSize = [length(x1d), length(y1d), length(z1d)];

Band = sub2ind(BandSearchSize, IJK(:,1), IJK(:,2), IJK(:,3));


% Create L, E, M
L = laplacian_3d_matrix(y1d,x1d,z1d, Lorder, Band);

Eplot = interp3_matrix(y1d, x1d, z1d, PointCloud.Location(:,2), PointCloud.Location(:,1), PointCloud.Location(:,3), porder, Band);
% Eplot = interpLagrange3D(BandSearchSize, MinPoint, PointCloud.Location, porder, Band, spacing);

Ecp = interp3_matrix(y1d, x1d, z1d, CP(:,2), CP(:,1), CP(:,3), porder, Band);
% Ecp = interpLagrange3D(BandSearchSize, MinPoint, CP, porder, Band, spacing);

M = lapsharp(L, Ecp);



% Extrapolate data to embedding

FaceInterpolateWeights = interpBarycenterTriangle(PointCloud, CP, CPFACE);

Signal = FaceInterpolateWeights *  PointCloud.Signal;



SignalImplicit = zeros(PointCloud.LocationCount, NumStepsImplicit);
SignalImplicit(:,1) = PointCloud.Signal;


WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsImplicit-1));



% G = make3DGaussianMatrix(y1d, x1d, z1d, sigma, spacing, Band, 4, 0);
% 
% GE =  G * Ecp ;
% 
% 
% 
% 
% for i = 1 : NumStepsImplicit - 1
%     
% %     SignalNew = Eplot * GE * FaceInterpolateWeights * Signal;
% %     
% %     SignalNew(abs(SignalNew)<10*eps) = 0;
% %     i
% %     SignalNew(SignalNew<0)
% %     
% %     SignalImplicit(:,i+1) = SignalNew;
%     SignalNew = GE * Signal;
%     
%     
%     SignalEplot = Eplot * SignalNew;
%     SignalEplot(abs(SignalEplot)<10*eps) = 0;
%     if ~isempty(SignalEplot(SignalEplot<0))
%         i
%         SignalEplot(SignalEplot<0)
%     end
%     
% %     SignalImplicit(:,i+1) = Eplot * SignalNew;
%     SignalImplicit(:,i+1) = SignalEplot;
%     
%     Signal = SignalNew;
%     
%     
%     
%     waitbar(i/NumStepsImplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsImplicit-1));
%     
% end


ItM = speye(size(M)) - tauImplicit * M;

for i = 1 : NumStepsImplicit - 1
    

    [SignalNew, flag, relres] = bicg(ItM, Signal, 1e-10, 30);

%     [SignalNew, flag, relres] = gmres(ItM, Signal, [], 1e-10, 30);

    if flag
        flag
        relres
    end
    
    
    Signal = SignalNew;
    

    
    SignalImplicit(:,i+1) = Eplot * SignalNew;
    
    waitbar(i/NumStepsImplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsImplicit-1));
    
end

waitbar(i/NumStepsImplicit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)



SignalImplicit(SignalImplicit<0) = 0;

figure
for i = 1 : NumStepsImplicit

    imshow(reshape(SignalImplicit(:,i),NumPointSurf,NumPointSurf),[]);
    drawnow
%     pause(0.01)
    
end

imshow(reshape(SignalImplicit(:,end),NumPointSurf,NumPointSurf),[]);


RadialDist3D = sqrt(sum(bsxfun(@minus, PointCloud.Location, PointCloud.Location(PointCloudCenter,:)).^2,2));

xgauss = 0:0.01:20;

for i = 1 : NumStepsImplicit
    clf
    plot(RadialDist3D, SignalImplicit(:,i),'r.')

    hold on
    axis([0 10 0 2*max(SignalImplicit(:,i))])
    
    
    t = (i-1)*sigma^2
    
    Gauss = exp(-(xgauss.^2) / (2*t));
    Gauss = Gauss * max(SignalImplicit(:,i));
    plot(xgauss, Gauss,'k:')


        
    pause

    
end





% Find the scale parameters
maxsample = 10; 
ws = 0 : 0.01 : maxsample;
NumSample = length(ws);
ScaleParameterSpatialImplicit = zeros(NumStepsImplicit,1);
cut = sqrt(log(2));
db3 = 1/sqrt(2);
ws2 = (ws.^2)';
H = ones(NumSample,1) ;
% h = 1 ./ ( ones(NumSample,1) + tauImplicit/spacing^2 * ws2 );
h = 1 ./ ( ones(NumSample,1) + tauImplicit * ws2 );

for i = 2 : NumStepsImplicit
    
    % Transfer function
    % 1st order
    H = H .* h;
    % Find the frequency at the cutoff values
    [uH, indH] = unique(H);
    CutoffFrequencyImplicit = interp1(uH, ws(indH), db3);
    % Change cutoff frequency to scale parameter
    ScaleParameterFrequencyImplicit = CutoffFrequencyImplicit / cut;   
%     ScaleParameterSpatialImplicit(i,1) =  spacing / ScaleParameterFrequencyImplicit;
    ScaleParameterSpatialImplicit(i,1) =  1 / ScaleParameterFrequencyImplicit;

end



% a = (hEdge/2) \ spacing

RadialDist3D = sqrt(sum(bsxfun(@minus, PointCloud.Location, PointCloud.Location(PointCloudCenter,:)).^2,2));

xgauss = 0:0.01:20;

for i = 1 : NumStepsImplicit
    clf
    plot(RadialDist3D, SignalImplicit(:,i),'r.')
    hold on
    axis([0 10 0 2*max(SignalImplicit(:,i))])
    
   
    Gauss = exp(-(xgauss.^2) / (2*ScaleParameterSpatialImplicit(i)^2));
    Gauss = Gauss * max(SignalImplicit(:,i));
    plot(xgauss, Gauss,'b-')
    
    
    i
    t = (i-1)*2*tauImplicit
    ScaleParameterSpatialImplicit(i)^2    
    
    
    Gauss = exp(-(xgauss.^2) / (2*t));
    Gauss = Gauss * max(SignalImplicit(:,i));
    plot(xgauss, Gauss,'k:')
   

    
    pause
    
end










%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures for paper %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% spacing = 0.1 [26 101 301]
% spacing = 0.25; [11 41 121]
% spacing = 0.5 [6 21 61]
for i = [5,15,35]
    figure
    imshow(reshape(SignalImplicit(:,i),NumPointSurf,NumPointSurf),[]);
    
    
end



RadialDist3D = sqrt(sum(bsxfun(@minus, PointCloud.Location, PointCloud.Location(PointCloudCenter,:)).^2,2));

xgauss = 0:0.01:20;



% spacing = 0.1 [26 101 301]
% spacing = 0.25; [11 41 121]
% spacing = 0.5 [6 21 61]
for i = [5,15,35]
    
    figure
    
    hold on
    
   
    Gauss = exp(-(xgauss.^2) / (2*ScaleParameterSpatialImplicit(i)^2));
    Gauss = Gauss * max(SignalImplicit(:,i));
    plot(xgauss, Gauss,'m--','linewidth',3)
    
    t = (i-1)*2*tauImplicit
    ScaleParameterSpatialImplicit(i)^2  
    
    
    Gauss = exp(-(xgauss.^2) / (2*t));
    Gauss = Gauss * max(SignalImplicit(:,i));
    plot(xgauss, Gauss,'b:','linewidth',3)
    
    plot(RadialDist3D, SignalImplicit(:,i),'r.','markersize',30)
    
    
    xlim([0 5])
    ylabel('Intensity')
    xlabel('Radial Distance')
    ax = gca;
    ax.XAxis.FontSize = 55;
    ax.YAxis.FontSize = 55;
    
%     pause
    
end
























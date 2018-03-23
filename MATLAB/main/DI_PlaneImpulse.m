% Andrew Rhodes
% ASEL
% March 2018


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
addpath(genpath('../'))
addpath('../src/')
addpath('../data/')
addpath(genpath('../models/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 2D image diffusion
MaxLevel = 70;
gausswindow = 2;
tau2D = 0.75;

MaxImageSize = 30;
MiddleImage = MaxImageSize / 2;


% 3D implicit surface (DI)
spacing = 0.25;
%     sigma = 0.25;
%     porder = 3;
Lorder = 2;
dim = 3;
tauImplicit = tau2D^2/2;
MaxTauImplicit = 5^2;
NumStepsImplicit = round(MaxTauImplicit / (2*tauImplicit) + 1);

alpha = 1;

MCporder = [2,3,4,5,6];

MCErrorDI = zeros(NumStepsImplicit,length(MCporder));

for MC = 1 : length(MCporder)
    
    porder = MCporder(MC);
    
    
    % % 2D image diffusion
    %     MaxLevel = 70;
    %     gausswindow = 2;
    %     tau2D = 0.75;
    %
    %     MaxImageSize = 30;
    %     MiddleImage = MaxImageSize / 2;
    %
    %
    %     % 3D implicit surface (DI)
    %     spacing = 0.25;
    %     sigma = 0.25;
    %     porder = 3;
    %     Lorder = 2;
    %     dim = 3;
    %     tauImplicit = tau2D^2/2;
    %     MaxTauImplicit = 5^2;
    %     NumStepsImplicit = round(MaxTauImplicit / (2*tauImplicit) + 1);
    
    
    
    % 3D explicit surface (DE)
    eSS = 1;
    % 3D explicit surface (DE)
    MaxSurfSize = 30;
    MaxTauExplicit  = 5^2;
    tauExplicit = tau2D^2/2;
    NumStepsExplcit = round(MaxTauExplicit / (2*tauExplicit) + 1);
    
    % If using Gaussian implicit surface diffusion, p=3 causes negative values
    % porder = 4
    
    NumPointSurf = length(1:eSS:MaxSurfSize+(1/eSS -1)*eSS);
    MiddleSurfPoint = round(NumPointSurf/2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup File Name Directions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    FileLocation = '../models/Plane/';
    FileName = strcat('Plane','_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.off');
    
    FileLocationCP = '../models/Plane/CPLaplace/';
    FileNameIJK = strcat('IJK','_p',num2str(porder),'_l',num2str(Lorder),'.mat');
    FileNameCP = strcat('CP','_p',num2str(porder),'_l',num2str(Lorder),'.mat');
    FileNameCPFACE = strcat('CPFACE','_p',num2str(porder),'_l',num2str(Lorder),'.mat');
    
    FileNameL = strcat('L','_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
    FileNameEplot = strcat('Eplot','_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
    FileNameEcp = strcat('Ecp','_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
    FileNameM = strcat('M','_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Explicit 3D Surface, Laplace-Beltrami, diffusion for impulse
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    if ~exist(fullfile(FileLocation, FileName), 'file')
        
        [PointCloud.Location, PointCloud.Face] = icosphere(NumberDivisions);
        
        save_off(PointCloud.Location, PointCloud.Face, fullfile(FileLocation, FileName))
        
        
    else
        
        [PointCloud.Location, PointCloud.Face] = read_off(fullfile(FileLocation, FileName));
        
        [m, n] = size(PointCloud.Location);
        if m < n
            PointCloud.Location = PointCloud.Location';
        end
        
        [m, n] = size(PointCloud.Face);
        if m < n
            PointCloud.Face = PointCloud.Face';
        end
        
    end
    
    PointCloudCenter = sub2ind([NumPointSurf, NumPointSurf], MiddleSurfPoint, MiddleSurfPoint);
    PointCloud.LocationCount = length(PointCloud.Location);
    PointCloud.FaceCount = length(PointCloud.Face);

    PointCloud.FaceArea = findFaceArea(PointCloud.Location, PointCloud.Face);
    PointCloud.Signal = zeros(PointCloud.LocationCount,1);
    PointCloud.Signal(PointCloudCenter,1) = 1;
    PointCloud = findMeshResolution(PointCloud, 'Model');
    
    
    
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
    
    
    
    if ~exist(fullfile(FileLocationCP,FileNameIJK), 'file')
        
        [IJK,DIST,CP,XYZ,CPFACE] = tri2cp(PointCloud.Face, PointCloud.Location, spacing, MinPoint, porder, Lorder/2);
        
        save(fullfile(FileLocationCP, FileNameIJK), 'IJK', '-v7.3')
        save(fullfile(FileLocationCP, FileNameCP), 'CP', '-v7.3')
        save(fullfile(FileLocationCP, FileNameCPFACE), 'CPFACE', '-v7.3')
        
    else
        
        load(fullfile(FileLocationCP, FileNameIJK))
        load(fullfile(FileLocationCP, FileNameCP))
        load(fullfile(FileLocationCP, FileNameCPFACE))
        
    end
      
    
    
    % XYZ = MinPoint - spacing + IJK * spacing
    
    BandSearchSize = [length(x1d), length(y1d), length(z1d)];
    
    Band = sub2ind(BandSearchSize, IJK(:,1), IJK(:,2), IJK(:,3));
    
    
    
    if ~exist(fullfile(FileLocationCP, FileNameL), 'file')
        
        % Create L, E, M
        L = laplacian_3d_matrix(y1d,x1d,z1d, Lorder, Band);
        
        Eplot = interp3_matrix(y1d, x1d, z1d, PointCloud.Location(:,2), PointCloud.Location(:,1), PointCloud.Location(:,3), porder, Band);
        % Eplot = interpLagrange3D(BandSearchSize, MinPoint, PointCloud.Location, porder, Band, spacing);
        
        Ecp = interp3_matrix(y1d, x1d, z1d, CP(:,2), CP(:,1), CP(:,3), porder, Band);
        % Ecp = interpLagrange3D(BandSearchSize, MinPoint, CP, porder, Band, spacing);
        
        M = lapsharp(L, Ecp);
        
        save(fullfile(FileLocationCP, FileNameL), 'L','-v7.3')
        save(fullfile(FileLocationCP, FileNameEplot), 'Eplot','-v7.3')
        save(fullfile(FileLocationCP, FileNameEcp), 'Ecp','-v7.3')
        save(fullfile(FileLocationCP, FileNameM), 'M','-v7.3')
        
    else
        load(fullfile(FileLocationCP, FileNameL))
        load(fullfile(FileLocationCP, FileNameEplot))
        load(fullfile(FileLocationCP, FileNameEcp))
        load(fullfile(FileLocationCP, FileNameM))
        
    end
    
    
    
    
    
    % Extrapolate data to embedding
    
    FaceInterpolateWeights = interpBarycenterTriangle(PointCloud, CP, CPFACE);
    
    Signal = FaceInterpolateWeights *  PointCloud.Signal;
    
    
    SignalImplicit = zeros(PointCloud.LocationCount, NumStepsImplicit);
    SignalImplicit(:,1) = PointCloud.Signal;
    
    
    WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsImplicit-1));
    
    ItM = speye(size(M)) - alpha * tauImplicit * M;
    
    for i = 1 : NumStepsImplicit - 1
        
        
        [SignalNew, flag, relres] = bicg(ItM, Signal, 1e-10, 60);
        
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
    
    %     figure
    %     for i = 1 : NumStepsImplicit
    %
    %         imshow(reshape(SignalImplicit(:,i),NumPointSurf,NumPointSurf),[]);
    %         drawnow
    %
    %     end
    
    
    
    ScaleParameterSpatialImplicit = findScaleParamter(tauImplicit, alpha, NumStepsImplicit, 1, 3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    RadialDist3D = sqrt(sum(bsxfun(@minus, PointCloud.Location, PointCloud.Location(PointCloudCenter,:)).^2,2));
    
    
    for i = 1 : NumStepsImplicit
        
        
        ErrorGauss = exp(-RadialDist3D.^2 / (2*ScaleParameterSpatialImplicit(i)^2));
        ErrorGauss = ErrorGauss * max(SignalImplicit(:,i));
        
        
        MCErrorDI(i,MC) = norm(ErrorGauss - SignalImplicit(:,i), inf);
        
        
    end
    
    
    
    
    
    
    
    
end




figure
loglog(1:NumStepsImplicit, MCErrorDI(:,1),'r--','linewidth',3)
hold on
loglog(1:NumStepsImplicit, MCErrorDI(:,2),'g:','linewidth',3)
loglog(1:NumStepsImplicit, MCErrorDI(:,3),'b-.','linewidth',3)
loglog(1:NumStepsImplicit, MCErrorDI(:,4),'m--','linewidth',3)
loglog(1:NumStepsImplicit, MCErrorDI(:,5),'k-','linewidth',3)








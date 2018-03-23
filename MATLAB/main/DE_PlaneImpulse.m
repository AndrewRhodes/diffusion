% Andrew Rhodes
% ASEL
% March 2018


close all
clear
clc

ProjectRoot = addprojectpaths % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 2D image diffusion
MaxLevel = 70;
gausswindow = 2;
tau2D = 0.75;

MaxImageSize = 30;
MiddleImage = MaxImageSize / 2;

alpha = 1;

MCeSS = [1, 0.5, 0.25, 0.1, 0.05];
% 3D explicit surface (DE)
MaxSurfSize = 30;
MaxTauExplicit  = 5^2;
tauExplicit = tau2D^2/2;
NumStepsExplicit = round(MaxTauExplicit / (2*tauExplicit) + 1);

MCErrorDE = zeros(NumStepsExplicit,length(MCeSS));

for MC = 1 : length(MCeSS)
    
    eSS = MCeSS(MC);
    
    % % 2D image diffusion
    % MaxLevel = 70;
    % gausswindow = 2;
    % tau2D = 0.75;
    %
    % MaxImageSize = 30;
    % MiddleImage = MaxImageSize / 2;
    
    
    % 3D implicit surface (DI)
    spacing = 0.25;
    sigma = 0.25;
    porder = 3;
    Lorder = 2;
    dim = 3;
    tauImplicit = tau2D^2/2;
    MaxTauImplicit = 5^2;
    NumStepsImplicit = round(MaxTauImplicit / (2*tauImplicit) + 1);
    
	ShowPlot = 1;
    
    % 3D explicit surface (DE)
    % eSS = 1; % Explicit Surface Spacing
    % MaxSurfSize = 30;
    % MaxTauExplicit  = 5^2;
    % tauExplicit = tau2D^2/2;
    % NumStepsExplcit = round(MaxTauExplicit / (2*tauExplicit) + 1);
    
    % If using Gaussian implicit surface diffusion, p=3 causes negative values
    % porder = 4
    
    
    NumPointSurf = length(1:eSS:MaxSurfSize+(1/eSS -1)*eSS);
    MiddleSurfPoint = round(NumPointSurf/2);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup File Name Directions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    FileLocation = '../models/Plane/';
    FileName = strcat('Plane','_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.off');
    
    FileLocationMeshLP = '../models/Plane/meshLP/';
    FileNameLapMat = strcat('LapMatMeshWeights','_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.mat');
    FileNameArea = strcat('Area','_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.mat');
    FileNamehEdge = strcat('hEdge2','_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.mat');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Explicit 3D Surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    if ~exist(fullfile(FileLocation, FileName), 'file')
        
        [xSurf3D, ySurf3D, zSurf3D] = ndgrid(1:eSS:MaxSurfSize+(1/eSS -1)*eSS,1:eSS:MaxSurfSize+(1/eSS -1)*eSS,0);
        
        PointCloud.Location = [xSurf3D(:), ySurf3D(:), zSurf3D(:)];
        PointCloud.Face = delaunay(xSurf3D(:), ySurf3D(:));
        
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
    % Create Laplace-Beltrami
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    if ~exist(fullfile(FileLocationMeshLP, FileNameLapMat), 'file')
        
        [LapMatMeshWeights, Area, hEdge2] = symmshlp_matrix(fullfile(FileLocation, FileName));
        
        save(fullfile(FileLocationMeshLP, FileNameLapMat), 'LapMatMeshWeights')
        save(fullfile(FileLocationMeshLP, FileNameArea), 'Area')
        save(fullfile(FileLocationMeshLP, FileNamehEdge), 'hEdge2')
    else
        
        load(fullfile(FileLocationMeshLP, FileNameLapMat))
        load( fullfile(FileLocationMeshLP, FileNameArea) )
        load( fullfile(FileLocationMeshLP, FileNamehEdge) )
        
    end
    
    % [LapMatMeshWeights, Area, hEdge2] = symmshlp_matrix('Plane120.off');
    hEdge = (hEdge2/2);
    
    A1 = sparse(1:length(Area),1:length(Area), 1./Area);
    
    LBM = A1 * LapMatMeshWeights;
    
    ItL = speye(length(Area),length(Area)) - alpha*tauExplicit * LBM;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Scale Paramter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ScaleParameterSpatialExplicit = findScaleParamter(tauExplicit, alpha, NumStepsExplicit, 1, 3);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform Diffusion and Error Analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	RadialDist3D = sqrt(sum(bsxfun(@minus, PointCloud.Location, PointCloud.Location(PointCloudCenter,:)).^2,2));
	GaussLocation = 0:0.01:max(PointCloud.Location)/2;
    if ShowPlot
        figure(1)
    end    
       
    SignalExplicit = zeros(PointCloud.LocationCount, NumStepsExplicit);
    SignalExplicit(:,1) = PointCloud.Signal;
    
    
    WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsExplicit-1));

    for i = 1 : NumStepsExplicit - 1
        
        [SignalExplicit(:,i+1), flag] = bicg(ItL, SignalExplicit(:,i), 1e-10, 60);
        if flag
            disp(flag)
        end
        
        waitbar(i/NumStepsExplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsExplicit-1));
        
        % Calculate Error
        TruthGauss = exp(-RadialDist3D.^2 / (2*ScaleParameterSpatialExplicit(i)^2));
        TruthGauss = TruthGauss * max(SignalExplicit(:,i));
        
        MCErrorDE(i,MC) = norm(TruthGauss - SignalExplicit(:,i), inf);
        
        % Show the plot
        if ShowPlot
            clf
            Gauss = exp(-GaussLocation.^2 / (2*ScaleParameterSpatialExplicit(i)^2));
            Gauss = Gauss * max(SignalExplicit(:,i));
            plot(GaussLocation, Gauss,'k')
            hold on
            plot(RadialDist3D, SignalExplicit(:,i),'r.')
            axis([0 10 0 2*max(SignalExplicit(:,i))])
        end
        
        
    end
    
    waitbar(i/NumStepsExplicit, WaitBar, sprintf('Diffusion Complete'));
    close(WaitBar)
    close(figure(1))
    
    
    
end

% save DE_Changing_eSS MCErrorDE

logx = [5,30];
logy = (10e-6).*logx.^(-2);

figure
loglog(1:NumStepsExplicit, MCErrorDE(:,1),'r--','linewidth',3)
hold on
loglog(1:NumStepsExplicit, MCErrorDE(:,2),'g:','linewidth',3)
loglog(1:NumStepsExplicit, MCErrorDE(:,3),'b-.','linewidth',3)
loglog(1:NumStepsExplicit, MCErrorDE(:,4),'m--','linewidth',3)
loglog(1:NumStepsExplicit, MCErrorDE(:,5),'k-','linewidth',3)
loglog(logx, logy,'k-','linewidth',3)

xl=get(gca,'XLim');
yl=get(gca,'YLim');
ht = text(6,10*10^(-8.5),'2nd Order');
set(ht,'Rotation',-11)
set(ht,'FontSize',50)


legend({'$\bar{\textbf{e}} = 1$','$\bar{\textbf{e}} = 0.5$','$\bar{\textbf{e}} = 0.25$','$\bar{\textbf{e}} = 0.1$','$\bar{\textbf{e}} = 0.05$'},'FontSize',40,'Interpreter','latex')
xlabel('Iterations')
ylabel('$\| error \|_{\infty}$','Interpreter','latex')
ax = gca;
ax.XAxis.FontSize = 55;
ax.YAxis.FontSize = 55;
xlim([10*10^(-1) 10e2])
ylim([10e-10 10e0])



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

MCeSS = [1, 0.5, 0.25, 0.1, 0.05]; % , 0.025, 0.01];
MCrho = [1, 2, 3, 4, 5, 6, 7, 8, 9];

MCError = zeros(length(MCeSS), 2, length(MCrho));
MCErrorAll = cell(length(MCeSS), length(MCrho));

for MCr = 1 : length(MCrho)

for MCe = 1 : length(MCeSS)
    
    clearvars -except MCr MCe MCeSS MCrho MCError MCErrorAll
    
    eSS = MCeSS(MCe)
    options.rho = MCrho(MCr)
    
%     options.dtype = 'geodesic';
    options.dtype = 'euclidean';
    
    ShowPlot = 0;
    
    % 2D image diffusion
%     MaxLevel = 70;
%     gausswindow = 2;
%     tau2D = 0.75;
    
%     MaxImageSize = 30;
%     MiddleImage = MaxImageSize / 2;
    


% tauImplicit = spacing / 8;
% MaxTauImplicit = 1/spacing;
% NumStepsImplicit = round(MaxTauImplicit);


    % 3D explicit surface (DE)
    tauExplicit = eSS / 8;
    MaxTauExplicit  = 1 / eSS;
    NumStepsExplicit = round(MaxTauExplicit) ;%round(MaxTauExplicit / (2*tauExplicit) + 1);
    if NumStepsExplicit == 1
        NumStepsExplicit = 2;
    end
	alpha = 1;
    
    
	MaxSurfSize = 30;
    NumPointSurf = length(1:eSS:MaxSurfSize+(1/eSS -1)*eSS);
    MiddleSurfPoint = round(NumPointSurf/2);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup File Name Directions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    FileLocation = '../models/Plane/';
    FileName = strcat('Plane','_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.off');
    
    FileLocationMeshLP = '../models/Plane/meshLP/';
    FileNameLapMat = strcat('LapMatMeshWeights','_eSS',num2str(eSS),'_',num2str(NumPointSurf),'_NumSigma',num2str(options.rho),'.mat');
    FileNameArea = strcat('Area','_eSS',num2str(eSS),'_',num2str(NumPointSurf),'_NumSigma',num2str(options.rho),'.mat');
    FileNamehEdge = strcat('hEdge2','_eSS',num2str(eSS),'_',num2str(NumPointSurf),'_NumSigma',num2str(options.rho),'.mat');
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Explicit 3D Surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%     
%     if ~exist(fullfile(FileLocation, FileName), 'file')
%         
%         [xSurf3D, ySurf3D, zSurf3D] = ndgrid(1:eSS:MaxSurfSize+(1/eSS -1)*eSS,1:eSS:MaxSurfSize+(1/eSS -1)*eSS,0);
%         
%         PointCloud.Location = [xSurf3D(:), ySurf3D(:), zSurf3D(:)];
%         PointCloud.Face = delaunay(xSurf3D(:), ySurf3D(:));
%         
%         save_off(PointCloud.Location, PointCloud.Face, fullfile(FileLocation, FileName))
% %         warning('Does Not Exist eSS = %0.3f', eSS)
%     else
% %          sprintf('Already Exist eSS = %0.3f',eSS)
%         [PointCloud.Location, PointCloud.Face] = read_off(fullfile(FileLocation, FileName));
%         
%         [m, n] = size(PointCloud.Location);
%         if m < n
%             PointCloud.Location = PointCloud.Location';
%         end
%         
%         [m, n] = size(PointCloud.Face);
%         if m < n
%             PointCloud.Face = PointCloud.Face';
%         end
%         
%     end
%     
% 
%         
%     PointCloudCenter = sub2ind([NumPointSurf, NumPointSurf], MiddleSurfPoint, MiddleSurfPoint);
%     
%     PointCloud.LocationCount = length(PointCloud.Location);
%     PointCloud.FaceCount = length(PointCloud.Face);
%     PointCloud.FaceArea = findFaceArea(PointCloud.Location, PointCloud.Face);
%     PointCloud.Signal = zeros(PointCloud.LocationCount,1);
%     PointCloud.Signal(PointCloudCenter,1) = 1;
%     PointCloud = findMeshResolution(PointCloud, 'Model');
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Laplace-Beltrami
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    if ~exist(fullfile(FileLocationMeshLP, FileNameLapMat), 'file')
        
%         [LapMatMeshWeights, Area, hEdge2] = symmshlp_matrix(fullfile(FileLocation, FileName), options);
%         
%         save(fullfile(FileLocationMeshLP, FileNameLapMat), 'LapMatMeshWeights','-v7.3')
%         save(fullfile(FileLocationMeshLP, FileNameArea), 'Area','-v7.3')
%         save(fullfile(FileLocationMeshLP, FileNamehEdge), 'hEdge2','-v7.3')
        
        warning('Does Not Exist eSS = %0.3f, Sigma = %i',eSS, options.rho)
    else
        
        sprintf('Already Exist eSS = %0.3f, Sigma = %i',eSS, options.rho)
        
%         load(fullfile(FileLocationMeshLP, FileNameLapMat))
%         load( fullfile(FileLocationMeshLP, FileNameArea) )
%         load( fullfile(FileLocationMeshLP, FileNamehEdge) )
        
    end
   
end
end
    
    hEdge = (hEdge2/2);
    
    AreaLength = length(Area);
    A1 = sparse(1:AreaLength,1:AreaLength, 1./Area);
    
    LBM = A1 * LapMatMeshWeights;
    
    ItL = speye(AreaLength,AreaLength) - alpha*tauExplicit * LBM;
    I23tL = speye(AreaLength,AreaLength) - (2/3)*alpha*tauExplicit * LBM;
    I611tL = speye(AreaLength,AreaLength) - (6/11)*alpha*tauExplicit * LBM;
    I1225tL = speye(AreaLength,AreaLength) - (12/25)*alpha*tauExplicit * LBM;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Scale Parameter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    ScaleParameterSpatialExplicit = findScaleParamter(tauExplicit, alpha, NumStepsExplicit, 1, 3);


    if max(ScaleParameterSpatialExplicit) > 5^2
        error('Too Much Diffusion, Decrease NumStepsExplicit')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Impulse diffusion on plane
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    RadialDist3D = sqrt(sum(bsxfun(@minus, PointCloud.Location, PointCloud.Location(PointCloudCenter,:)).^2,2));
    GaussLocation = 0:0.01:max(PointCloud.Location)/2;
    if ShowPlot
        figure(1)
    end
    
    
    SignalExplicit = zeros(PointCloud.LocationCount, NumStepsExplicit);
    SignalExplicit(:,1) = PointCloud.Signal;
    
    AbsErr = zeros(NumStepsExplicit,1);
   
	WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsExplicit-1));
    for i = 1 : NumStepsExplicit - 1
        
        % Perform Diffusion
        if i == 1
% % %         [SignalExplicit(:,i+1), flag] = bicg(ItL, SignalExplicit(:,i), 1e-10, 60);
        [SignalExplicit(:,i+1), flag] = gmres(ItL, SignalExplicit(:,i), [], 1e-10, 100);
        elseif i == 2
            [SignalExplicit(:,i+1), flag] = gmres(I23tL, (4/3)*SignalExplicit(:,i) - (1/3)*SignalExplicit(:,i-1), [], 1e-10, 100);
        elseif i == 3
            [SignalExplicit(:,i+1), flag] = gmres(I611tL, (18/11)*SignalExplicit(:,i) - (9/11)*SignalExplicit(:,i-1) + (2/11)*SignalExplicit(:,i-2), [], 1e-10, 100);
        else
            [SignalExplicit(:,i+1), flag] = gmres(I1225tL, (48/25)*SignalExplicit(:,i)-(36/25)*SignalExplicit(:,i-1) + (16/25)*SignalExplicit(:,i-2) - (3/25)*SignalExplicit(:,i-3), [], 1e-10, 100);
        end
       
        
        if flag
            disp(flag)
        end
        
        waitbar(i/NumStepsExplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsExplicit-1));
        
        % Calculate Error
        TruthGauss = exp(-RadialDist3D.^2 / (2*ScaleParameterSpatialExplicit(i)^2));
        TruthGauss = TruthGauss * max(SignalExplicit(:,i));
        
        AbsErr(i+1,1) = norm(TruthGauss - SignalExplicit(:,i), inf);
        
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
        
        pause(0.1)
    end
    
    waitbar(i/NumStepsExplicit, WaitBar, sprintf('Diffusion Complete'));
    close(WaitBar)
    close(figure(1))
      
    
    MCError(MCe, 1:2, MCr) = [NumStepsExplicit, AbsErr(NumStepsExplicit - 1)];
    MCErrorAll{MCe, MCr} = [(1:NumStepsExplicit)', AbsErr];
    
    
end
end


colors = ['c','m','k','b','r','c','m','g','k'];
shapes = {'','','.','^','.','*','+','o','s'};
lin  = {'-','-',':','-','-','-','-.',':','--'};
msizes = [8, 8, 8, 8, 8, 8, 10, 12, 15];


Points = zeros(size(MCErrorAll,1),2, size(MCErrorAll,2));
for j = 1 : size(MCErrorAll,2)
    for i = 1 : size(MCErrorAll,1)
        Points(i,1:2,j) = [MCErrorAll{i,j}(end,1), MCErrorAll{i,j}(end,2)];
    end
end


figure('units','normalized','outerposition',[0 0 1 1])
for i = 1 : size(Points,3)
    loglog( Points(:,1,i), Points(:,2,i), strcat(colors(i), shapes{i}, lin{i}), 'linewidth', 3, 'markersize', msizes(i) )
    hold on
%     pause
end
logx = [1,100];
logy = (10e1).*logx.^(-6);
loglog(logx, logy,'k-')




% logx = [5,30];
% logy = (10e-6).*logx.^(-2);
% 
% figure
% loglog(1:NumStepsExplicit, MCErrorDE(:,1),'r--','linewidth',3)
% hold on
% loglog(1:NumStepsExplicit, MCErrorDE(:,2),'g:','linewidth',3)
% loglog(1:NumStepsExplicit, MCErrorDE(:,3),'b-.','linewidth',3)
% loglog(1:NumStepsExplicit, MCErrorDE(:,4),'m--','linewidth',3)
% loglog(1:NumStepsExplicit, MCErrorDE(:,5),'k-','linewidth',3)
% loglog(logx, logy,'k-','linewidth',3)
% 
% xl=get(gca,'XLim');
% yl=get(gca,'YLim');
% ht = text(6,10*10^(-8.5),'2nd Order');
% set(ht,'Rotation',-11)
% set(ht,'FontSize',50)
% 
% 
% legend({'$\bar{\textbf{e}} = 1$','$\bar{\textbf{e}} = 0.5$','$\bar{\textbf{e}} = 0.25$','$\bar{\textbf{e}} = 0.1$','$\bar{\textbf{e}} = 0.05$'},'FontSize',40,'Interpreter','latex')
% xlabel('Iterations')
% ylabel('$\| error \|_{\infty}$','Interpreter','latex')
% ax = gca;
% ax.XAxis.FontSize = 55;
% ax.YAxis.FontSize = 55;
% xlim([10*10^(-1) 10e2])
% ylim([10e-10 10e0])



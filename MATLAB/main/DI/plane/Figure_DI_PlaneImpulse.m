% Andrew Rhodes
% ASEL
% March 2018


close all
clear
clc

addpath('../../src/')
ProjectRoot = setupprojectpaths % Additional Paths


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCporder = [2, 3, 4, 5, 6, 7];
MCspacing = [1, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01];
%
MCError = zeros(length(MCspacing), 2, length(MCporder));
MCErrorAll = cell(length(MCspacing), length(MCporder));

for MCp = 1 : length(MCporder)
    for MCs = 1 : length(MCspacing)
        
        clearvars -except MCp MCs MCporder MCspacing MCError MCErrorAll ProjectRoot
        warning on
        
        porder = MCporder(MCp)
        spacing = MCspacing(MCs)
        
        ShowPlot = 0;
        
        tau2D = 0.75;
        eSS = 1;
        % 3D explicit surface (DE)
        MaxSurfSize = 30;
        NumPointSurf = length(1:eSS:MaxSurfSize+(1/eSS -1)*eSS);
        MiddleSurfPoint = round(NumPointSurf/2);
        
        
        
        % 3D implicit surface (DI)
        %     spacing = 0.25;
        Lorder = 2;
        dim = 3;
        tauImplicit = (spacing/8);
        MaxTauImplicit = 10/spacing;
        NumStepsImplicit = round(MaxTauImplicit);
        alpha = 1;
        
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup File Name Directions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FileLocation = strcat(ProjectRoot,'/models/Plane/');
        FileName = strcat('Plane','_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.off');
        
        FileLocationCP = strcat(ProjectRoot,'/models/Plane/CPLaplace/');
        FileNameIJK = strcat('IJK','_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.mat');
        FileNameCP = strcat('CP','_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.mat');
        FileNameCPFACE = strcat('CPFACE','_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.mat');
        
        FileNameL = strcat('L','_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.mat');
        FileNameEplot = strcat('Eplot','_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.mat');
        FileNameEcp = strcat('Ecp','_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.mat');
        FileNameM = strcat('M','_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'_eSS',num2str(eSS),'_',num2str(NumPointSurf),'.mat');
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create Explicit 3D Surface, Laplace-Beltrami, diffusion for impulse
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        if ~exist(fullfile(FileLocation, FileName), 'file')
            warning('Does Not Exist eSS = %0.3f, NumPointSurf = %i',eSS, NumPointSurf)
            
            [xSurf3D, ySurf3D, zSurf3D] = ndgrid(1:eSS:MaxSurfSize+(1/eSS -1)*eSS,1:eSS:MaxSurfSize+(1/eSS -1)*eSS,0);
            
            PointCloud.Location = [xSurf3D(:), ySurf3D(:), zSurf3D(:)];
            PointCloud.Face = delaunay(xSurf3D(:), ySurf3D(:));
            
            save_off(PointCloud.Location, PointCloud.Face, fullfile(FileLocation, FileName))
            
        else
            sprintf('Already Exist eSS = %0.3f, NumPointSurf = %i',eSS, NumPointSurf)
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
        if porder < 2
            bandwidth = 1.00001*spacing*sqrt((dim-1)*((2+1)/2)^2 + ((Lorder/2+(2+1)/2)^2));
        else
            bandwidth = 1.00001*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));
        end
        % Implicit Clostest Point Method
        MinPoint = round(min(PointCloud.Location) - bandwidth - spacing, 1);
        MaxPoint = round(max(PointCloud.Location) + bandwidth + spacing, 1);
        
        % Gaussian Method
        % MinPoint = round(min(PointCloud.Location) - bandwidth - 2*spacing, 1);
        % MaxPoint = round(max(PointCloud.Location) + bandwidth + 2*spacing, 1);
        
        x1d = (MinPoint(1):spacing:MaxPoint(1))';
        y1d = (MinPoint(2):spacing:MaxPoint(2))';
        z1d = (MinPoint(3):spacing:MaxPoint(3))';
        
        
        
        if ~exist(fullfile(FileLocationCP,FileNameIJK), 'file') || ~exist(fullfile(FileLocationCP,FileNameCP), 'file') || ~exist(fullfile(FileLocationCP,FileNameCPFACE), 'file')
%             warning('Does Not Exist IJK, CP, CPFACE, Lambda = %0.3f, Porder = %i', spacing, porder)
            
            [IJK,DIST,CP,XYZ,CPFACE] = tri2cp(PointCloud.Face, PointCloud.Location, spacing, MinPoint, porder, Lorder/2);
            
            save(fullfile(FileLocationCP, FileNameIJK), 'IJK', '-v7.3')
            save(fullfile(FileLocationCP, FileNameCP), 'CP', '-v7.3')
            save(fullfile(FileLocationCP, FileNameCPFACE), 'CPFACE', '-v7.3')
            
        else
%             sprintf('Already Exist IJK, CP, CPFACE, Lambda = %0.3f, Porder = %i',spacing, porder)
            
            load(fullfile(FileLocationCP, FileNameIJK))
            load(fullfile(FileLocationCP, FileNameCP))
            load(fullfile(FileLocationCP, FileNameCPFACE))
            
        end
        
        
        
        % XYZ = MinPoint - spacing + IJK * spacing
        
        BandSearchSize = [length(x1d), length(y1d), length(z1d)];
        
        Band = sub2ind(BandSearchSize, IJK(:,1), IJK(:,2), IJK(:,3));
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create Laplace-Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if ~exist(fullfile(FileLocationCP, FileNameL), 'file') || ~exist(fullfile(FileLocationCP, FileNameEplot), 'file') || ~exist(fullfile(FileLocationCP, FileNameEcp), 'file') || ~exist(fullfile(FileLocationCP, FileNameM), 'file')
            warning('Does Not Exist L,Ecp,M Lambda = %0.3f, Porder = %i',spacing, porder)

            % Create L, E, M
            L = laplacian_3d_matrix(x1d, y1d, z1d, Lorder, Band);
            
            Eplot = interp3_matrix(x1d, y1d, z1d, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), porder, Band);
            % Eplot = interpLagrange3D(BandSearchSize, MinPoint, PointCloud.Location, porder, Band, spacing);
                        
            Ecp = interp3_matrix(x1d, y1d, z1d, CP(:,1), CP(:,2), CP(:,3), porder, Band);
            % Ecp = interpLagrange3D(BandSearchSize, MinPoint, CP, porder, Band, spacing);
            

            M = lapsharp(L, Ecp);
            
            save(fullfile(FileLocationCP, FileNameL), 'L','-v7.3')
            save(fullfile(FileLocationCP, FileNameEplot), 'Eplot','-v7.3')
            save(fullfile(FileLocationCP, FileNameEcp), 'Ecp','-v7.3')
            save(fullfile(FileLocationCP, FileNameM), 'M','-v7.3')
            
        else
            
            sprintf('Already Exist L,Ecp,M Lambda = %0.3f, Porder = %i',spacing, porder)
                    load(fullfile(FileLocationCP, FileNameL))
                    load(fullfile(FileLocationCP, FileNameEplot))
                    load(fullfile(FileLocationCP, FileNameEcp))
                    load(fullfile(FileLocationCP, FileNameM))
            
        end
        

    
    ItM = speye(size(M)) - alpha * tauImplicit * M;
    
%     I23tM = speye(size(M)) - (2/3)*alpha*tauImplicit * M;
% 
%     I611M = speye(size(M)) - (6/11)*alpha*tauImplicit * M;
% 
%     I1225M = speye(size(M)) - (12/25)*alpha*tauImplicit * M;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Scale Parameter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ScaleParameterSpatialImplicit = findScaleParamter(tauImplicit, alpha, NumStepsImplicit, 1, 3);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Performe Diffusion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    RadialDist3D = sqrt(sum(bsxfun(@minus, PointCloud.Location, PointCloud.Location(PointCloudCenter,:)).^2,2));
    GaussLocation = 0:0.01:max(PointCloud.Location)/2;
    
    
    % Extrapolate data to embedding
    
    FaceInterpolateWeights = interpBarycenterTriangle(PointCloud, CP, CPFACE);
    
    Signal = FaceInterpolateWeights *  PointCloud.Signal;
    
    SignalImplicit = zeros(PointCloud.LocationCount, NumStepsImplicit);
    SignalImplicit(:,1) = PointCloud.Signal;
    
        
    NumIter = min(length(Band),100);
    
    if ShowPlot
        figure(1)
    end
    
	AbsErr = zeros(NumStepsImplicit,1);


	WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsImplicit-1));
    for i = 1 : NumStepsImplicit - 1
        
        %     if i == 1
%         Signal(:,i+1) = ItM \ Signal(:,i);
%         [SignalNew, flag, relres] = bicg(ItM, Signal, 1e-10, 60);
        [SignalNew, flag] = gmres(ItM, Signal, [], 1e-10, NumIter);
%     elseif i == 2
% %         Signal(:,i+1) = I23tM \ ((4/3)*Signal(:,i) - (1/3)*Signal(:,i-1));
%         [Signal(:,i+1), flag] = gmres(I23tM, (4/3)*Signal(:,i) - (1/3)*Signal(:,i-1), [], 1e-10, NumIter);
%     elseif i == 3
% %         Signal(:,i+1) = I611M \ ((4/3)*Signal(:,i) - (1/3)*Signal(:,i-1));
%         [Signal(:,i+1), flag] = gmres(I611M, (18/11)*Signal(:,i) - (9/11)*Signal(:,i-1) + (2/11)*Signal(:,i-2), [], 1e-10, NumIter);    
%     else
% %         Signal(:,i+1) = I1225M \ ((4/3)*Signal(:,i) - (1/3)*Signal(:,i-1));
%         [Signal(:,i+1), flag] = gmres(I1225M, (48/25)*Signal(:,i)-(36/25)*Signal(:,i-1) + (16/25)*Signal(:,i-2) - (3/25)*Signal(:,i-3), [], 1e-10, NumIter);    
%     end
        
               
        if flag
            disp(flag)         
        end
        
        % Rese the signal
        Signal = SignalNew;
        SignalImplicit(:,i+1) = Eplot * SignalNew;
        SignalImplicit( SignalImplicit(:,i+1) < 0, i+1) = 0;

        % Calculat Error
        TruthGauss = exp(-RadialDist3D.^2 / (2*ScaleParameterSpatialImplicit(i)^2));
        TruthGauss = TruthGauss * max(SignalImplicit(:,i));
        
        AbsErr(i+1,1) = norm(TruthGauss - SignalImplicit(:,i), inf);
        
        % Show the plot
        if ShowPlot
%             imshow(reshape(SignalImplicit(:,i),NumPointSurf,NumPointSurf),[]);
%             drawnow
%             pause(0.1)

            
            clf
            Gauss = exp(-GaussLocation.^2 / (2*ScaleParameterSpatialImplicit(i)^2));
            Gauss = Gauss * max(SignalImplicit(:,i));
            plot(GaussLocation, Gauss,'k')
            hold on
            plot(RadialDist3D, SignalImplicit(:,i),'r.')
            axis([0 10 0 2*max(SignalImplicit(:,i))])
            pause
            
        end
        
        waitbar(i/NumStepsImplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsImplicit-1));
        
    end
    
    waitbar(i/NumStepsImplicit, WaitBar, sprintf('Diffusion Complete'));
    close(WaitBar)
    close(figure(1))
    
    MCError(MCs, 1:2, MCp) = [NumStepsImplicit, AbsErr(NumStepsImplicit - 1)];
    MCErrorAll{MCs, MCp} = [(1:NumStepsImplicit)', AbsErr];
    
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
logx = [10,50];
logy = (10*10^(-0.2)).*logx.^(-3);
loglog(logx, logy,'k-')


% figure
% loglog(1:NumStepsImplicit, MCErrorDI(:,1),'r--','linewidth',3)
% hold on
% loglog(1:NumStepsImplicit, MCErrorDI(:,2),'g:','linewidth',3)
% loglog(1:NumStepsImplicit, MCErrorDI(:,3),'b-.','linewidth',3)
% loglog(1:NumStepsImplicit, MCErrorDI(:,4),'m--','linewidth',3)
% loglog(1:NumStepsImplicit, MCErrorDI(:,5),'k-','linewidth',3)








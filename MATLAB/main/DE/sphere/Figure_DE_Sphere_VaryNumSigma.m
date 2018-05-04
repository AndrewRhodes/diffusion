% Andrew Rhodes
% ASEL
% March 2018


% close all
clear
clc

addpath('../../src/')
ProjectRoot = setupprojectpaths % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCdiv = [1, 2, 3, 4, 5, 6];
MCrho = [ 3, 4, 5, 6, 7];

MCError = zeros(length(MCdiv), 2, length(MCrho));
MCErrorAll = cell(length(MCdiv), length(MCrho));

for MCr = 1 : length(MCrho)
    
    for MCd = 1 : length(MCdiv)
        
        clearvars -except MCr MCd MCdiv MCrho MCErrorAll MCError ProjectRoot
        
        NumberDivisions = MCdiv(MCd)
        options.rho = MCrho(MCr)
        
        
        options.dtype = 'geodesic';
%         options.dtype = 'euclidean';
        
        alpha = 1;
        
        ShowPlot = 1;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup File Name Directions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FileLocation = strcat(ProjectRoot,'/models/sphere/');
        FileName = strcat('Icosphere',num2str(NumberDivisions),'.off');
        
        FileLocationMeshLP = strcat(ProjectRoot,'/models/sphere/meshLP/');
        FileNameLapMat = strcat('LapMatMeshWeights','_Div',num2str(NumberDivisions),'_NumSigma',num2str(options.rho),'_',options.dtype,'.mat');
        FileNameArea = strcat('Area','_Div',num2str(NumberDivisions),'_NumSigma',num2str(options.rho),'_',options.dtype,'.mat');
        FileNamehEdge = strcat('hEdge2','_Div',num2str(NumberDivisions),'_NumSigma',num2str(options.rho),'_',options.dtype,'.mat');
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make the Sphere
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if ~exist(fullfile(FileLocation, FileName), 'file')
            
            [Location, Faces] = icosphere(NumberDivisions);
            [VerticesOut, FacesOut] = clearMeshDuplicates(Location, Faces );
            
            save_off(VerticesOut, FacesOut, fullfile(FileLocation, FileName))
            
            Sphere.Face = FacesOut;
            Sphere.Location = VerticesOut;
            
        else
            
            [Sphere.Location, Sphere.Face] = read_off(fullfile(FileLocation, FileName));
            
            [m, n] = size(Sphere.Location);
            if m < n
                Sphere.Location = Sphere.Location';
            end
            
            [m, n] = size(Sphere.Face);
            if m < n
                Sphere.Face = Sphere.Face';
            end
            
        end
        
        
        
        Sphere.FaceCount = size(Sphere.Face, 1);
        Sphere.LocationCount = size(Sphere.Location,1);
        Sphere.FaceArea = findFaceArea(Sphere.Location,Sphere.Face);
        
        Sphere.FaceArea = findFaceArea(Sphere.Location, Sphere.Face);
        Sphere.Signal = zeros(Sphere.LocationCount,1);
        Sphere = findMeshResolution(Sphere, 'Model');
        
        
        % % % % % % % % % %
        tauExplicit = Sphere.Resolution / 8;
        MaxTauExplicit = 1 / Sphere.Resolution;
        NumStepsExplicit = round(MaxTauExplicit);
        % % % % % % % % % %
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make the Signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [Sphere.Theta, Sphere.Phi, Sphere.Radius] = cart2sph(Sphere.Location(:,1) ,Sphere.Location(:,2), Sphere.Location(:,3));
        
        ExactSignal = @(sigma, Phi) exp(-sigma^2)*sin(Phi);
        
        SignalOriginal = ExactSignal(0, Sphere.Phi);
        Sphere.Signal = SignalOriginal;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Construct the Laplace Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        if ~exist(fullfile(FileLocationMeshLP, FileNameLapMat), 'file') || ~exist(fullfile(FileLocationMeshLP, FileNameArea), 'file') || ~exist(fullfile(FileLocationMeshLP, FileNamehEdge), 'file')
            
            [LapMatMeshWeights, Area, hEdge2] = symmshlp_matrix(fullfile(FileLocation, FileName), options);
            
            save(fullfile(FileLocationMeshLP, FileNameLapMat), 'LapMatMeshWeights', '-v7.3')
            save(fullfile(FileLocationMeshLP, FileNameArea), 'Area', '-v7.3')
            save(fullfile(FileLocationMeshLP, FileNamehEdge), 'hEdge2', '-v7.3')
        else
            
            load( fullfile(FileLocationMeshLP, FileNameLapMat))
            load( fullfile(FileLocationMeshLP, FileNameArea) )
            load( fullfile(FileLocationMeshLP, FileNamehEdge) )
            
        end
        
        hEdge = (hEdge2/2);
                
        Alength = length(Area);
        
        A1 = sparse(1:Alength, 1:Alength, 1./Area);
        
        LBM = A1 * LapMatMeshWeights;
        
        ItL = speye(Alength, Alength) - alpha*tauExplicit * LBM;
        
        I23tL = speye(Alength, Alength) - (2/3)*alpha*tauExplicit * LBM;
        
        I611tL = speye(Alength, Alength) - (6/11)*alpha*tauExplicit * LBM;
        
        I1225tL = speye(Alength, Alength) - (12/25)*alpha*tauExplicit * LBM;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Scale Parameter Estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        ScaleParameter = findScaleParamter(tauExplicit, alpha, NumStepsExplicit, 'natural', '3d');
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform Diffusion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        SignalExplicit = zeros(Sphere.LocationCount, NumStepsExplicit);
        SignalExplicit(:,1) = Sphere.Signal;
        
        AbsErr = zeros(NumStepsExplicit,1);
        
        
        WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsExplicit-1));
        
        if ShowPlot
            figure(1)
        end
        for i = 1 : NumStepsExplicit - 1
            
            
                if i == 1
                    [SignalExplicit(:,i+1), flag] = bicg(ItL, SignalExplicit(:,i), 1e-10, 100);
                elseif i == 2
                    [SignalExplicit(:,i+1), flag] = bicg(I23tL, (4/3)*SignalExplicit(:,i) - (1/3)*SignalExplicit(:,i-1), 1e-10, 100);
                elseif i == 3
                    [SignalExplicit(:,i+1), flag] = bicg(I611tL, (18/11)*SignalExplicit(:,i) - (9/11)*SignalExplicit(:,i-1) + (2/11)*SignalExplicit(:,i-2), 1e-10, 100);
                else
                    [SignalExplicit(:,i+1), flag] = bicg(I1225tL, (48/25)*SignalExplicit(:,i) - (36/25)*SignalExplicit(:,i-1) + (16/25)*SignalExplicit(:,i-2) - (3/25)*SignalExplicit(:,i-3), 1e-10, 100);
                end
            
            
            if flag
                disp(flag)
            end
            
            Truth = ExactSignal(ScaleParameter(i+1), Sphere.Phi);
            
            AbsErr(i+1,1) = norm(Truth - SignalExplicit(:,i+1), inf);
            
            if ShowPlot
                clf
                plot(Sphere.Phi, SignalOriginal,'ko')
                hold on
                plot(Sphere.Phi, Truth, 'gd')
                plot(Sphere.Phi, SignalExplicit(:,i+1),'r.')
            end
            
            
            waitbar(i/NumStepsExplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsExplicit-1));
        end
        
        if ShowPlot
            close(figure(1))
        end
        
        waitbar(i/NumStepsExplicit, WaitBar, sprintf('Diffusion Complete'));
        close(WaitBar)
        
        
        
        MCError(MCd, 1:2, MCr) = [NumStepsExplicit, AbsErr(NumStepsExplicit - 1)];
        MCErrorAll{MCd, MCr} = [(1:NumStepsExplicit)', AbsErr];
        
        
    end
    
end




colors = ['k','k','b','m','k','c','g','r'];
shapes = {'','','^','p','s','o','*','+'};
lin  = {':','-.','--','-','-','-','-.',':','--'};
msizes = [8, 8, 8, 8, 17, 14, 12, 8];

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
end




% % % 1st Order Line
logx = [100,500];
logy = (10*10^(0)).*logx.^(-2);
loglog(logx, logy,'k-','linewidth',3)

% xl=get(gca,'XLim');
% yl=get(gca,'YLim');
text1 = text(15,10*10^(-4.4),'$1^{st}$ Order','Interpreter','latex');
set(text1,'Rotation',-16)
set(text1,'FontSize',50)



% % % % 0th Order Line
% logx = [10,50];
% logy = (10*10^(-1.5)).*logx.^(0);
% loglog(logx, logy,'k-','linewidth',3)
%
% % xl=get(gca,'XLim');
% % yl=get(gca,'YLim');
% text1 = text(15,10*10^(-1.4),'$0^{th}$ Order','Interpreter','latex');
% set(text1,'Rotation',0)
% set(text1,'FontSize',50)

% % % Axis Labels
xlabel('N')
ylabel('$\| $error$ \|_{\infty}$','Interpreter','latex')
ax = gca;
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;



xlim([10*10^(-0.75) 10*10^(1.4)])

xticks([10e0 10e1 10e2 10e3])
yticks([10e-10, 10e-8, 10e-6, 10e-5,10e-4,10e-3, 10e-2, 10e0, 10e2, 10e4])
yticklabels({'10e-10','10e-8','10e-6','10e-5','10e-4','10e-3', '10e-2', '10e0', '10e2', '10e4'})


hleg = legend({'$m_\sigma = 1$','$m_\sigma = 2$','$m_\sigma = 3$','$m_\sigma = 4$','$m_\sigma = 5$','$m_\sigma = 6$','$m_\sigma = 7$','$m_\sigma = 8$'},'Interpreter','latex');
set(hleg, 'position', [0.80 0.58 0.07 0.2], 'FontSize', 35)












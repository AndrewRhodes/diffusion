% Andrew Rhodes
% ASEL
% March 2018

% Diffusion accuracy comparison for a signal on a sphere.


close all
clear
clc

addpath('../../src/')
ProjectRoot = setupprojectpaths; % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCspacing = [0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001];

MCnumsigma = [1, 2, 3, 4, 5, 6, 7, 8, 9];

MCError = zeros(length(MCspacing), 2, length(MCnumsigma));
MCErrorAll = cell(length(MCspacing), length(MCnumsigma));


for MCsigma = 1 : length(MCnumsigma)
    
    for MCs = 1 : length(MCspacing)
        
        
        clearvars -except MCspacing MCs MCsigma MCError MCErrorAll MCnumsigma ProjectRoot
        spacing = MCspacing(MCs)
        numsigmas = MCnumsigma(MCsigma)
        
        alpha = 1;
        porder = 4;
        dim = 2;
%         Lorder = 2;
        % spacing = 0.01;
        % sigma <= spacing
        sigma = spacing;
        % numsigmas = 7;
        LimitFarPoints = 1;
        
        if spacing > sigma
            bandwidth = 1.002*spacing*sqrt((dim-1)*((porder+1)/2)^2 + (((numsigmas*sigma)/spacing+(porder+1)/2)^2));
        else
            bandwidth = 1.002*numsigmas*sigma*sqrt((dim-1)*((porder+1)/2)^2 + (((numsigmas*sigma)/spacing+(porder+1)/2)^2));
        end
        ShowPlot = 0;
        NumCirclePoints = 1000;
        
        % tauImplicit = spacing / 8;
        % MaxTauImplicit = 1/spacing;
        NumSteps = round(1/spacing); % 5000 ; %
        
        ExactSignal = @(sigma, theta) exp(-sigma^2/2)*cos(theta) + exp(-9*sigma^2/2)*sin(3*theta);
%         ExactSignal = @(sigma, theta) exp(-sigma^2/2)*cos(theta);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make the Circle and Signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Theta = linspace(0, 2*pi, NumCirclePoints)';
        
        Radius = ones(size(Theta));
        [xp,yp] = pol2cart(Theta, Radius);
        Circle.Location(:,1) = xp(:);
        Circle.Location(:,2) = yp(:);
        Circle.LocationCount = length(Circle.Location);
        
        
        SignalOriginal = ExactSignal(0, Theta);
        
        Circle.Signal = SignalOriginal;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup File Name Directions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FileLocationCPGauss = strcat(ProjectRoot,'/models/Circle/CPGauss/');
        FileNameCP = strcat('CP','_NumPoints',num2str(NumCirclePoints),'_p',num2str(porder),'_s',num2str(spacing),'_NumSigma',num2str(numsigmas),'.mat');
        FileNameDIST = strcat('DIST','_NumPoints',num2str(NumCirclePoints),'_p',num2str(porder),'_s',num2str(spacing),'_NumSigma',num2str(numsigmas),'.mat');
        
        FileNameEcp = strcat('Ecp','_NumPoints',num2str(NumCirclePoints),'_p',num2str(porder),'_s',num2str(spacing),'_NumSigma',num2str(numsigmas),'.mat');
        FileNameEplot = strcat('Eplot','_NumPoints',num2str(NumCirclePoints),'_p',num2str(porder),'_s',num2str(spacing),'_NumSigma',num2str(numsigmas),'.mat');
        FileNameGCart = strcat('GCart','_NumPoints',num2str(NumCirclePoints),'_p',num2str(porder),'_s',num2str(spacing),'_NumSigma',num2str(numsigmas),'_LimitFar',num2str(LimitFarPoints),'.mat');
        FileNameG = strcat('G','_NumPoints',num2str(NumCirclePoints),'_p',num2str(porder),'_s',num2str(spacing),'_NumSigma',num2str(numsigmas),'_LimitFar',num2str(LimitFarPoints),'.mat');
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Construct the Laplace Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        MinPoint = min(Circle.Location) - bandwidth - ceil(numsigmas * sigma);
        MaxPoint = max(Circle.Location) + bandwidth + ceil(numsigmas * sigma);
        
        % MinPoint = round(min(Circle.Location(:,1)) - bandwidth, 1);
        % MaxPoint = round(max(Circle.Location(:,1)) + bandwidth, 1);
        
        
        x1d = (MinPoint(1):spacing:MaxPoint(1))';
        y1d = (MinPoint(2):spacing:MaxPoint(2))';
        
        [GridX, GridY] = meshgrid(x1d, y1d);
        
        if ~exist( fullfile(FileLocationCPGauss, FileNameCP), 'file') || ~exist( fullfile(FileLocationCPGauss, FileNameDIST), 'file')
            [CP(:,1), CP(:,2), DIST] = cpCircle(GridX(:), GridY(:));
           
            save(fullfile(FileLocationCPGauss, FileNameCP), 'CP', '-v7.3')
            save(fullfile(FileLocationCPGauss, FileNameDIST), 'DIST', '-v7.3')
        else
            load(fullfile(FileLocationCPGauss, FileNameCP), 'CP')
            load(fullfile(FileLocationCPGauss, FileNameDIST), 'DIST')
        end
        
        
        band = find(abs(DIST) <= bandwidth);
        
        CP = CP(band,:);
        
        GridXBand = GridX(band);
        GridYBand = GridY(band);
        
        
        [CPTheta, CPr] = cart2pol(GridXBand,GridYBand);
        
        
        CPSignal = ExactSignal(0, CPTheta);
        
        
        if ~exist( fullfile(FileLocationCPGauss, FileNameEcp), 'file') || ~exist( fullfile(FileLocationCPGauss, FileNameEplot), 'file') || ~exist( fullfile(FileLocationCPGauss, FileNameGCart), 'file') || ~exist( fullfile(FileLocationCPGauss, FileNameG), 'file')

            Ecp = interp2_matrix(x1d, y1d, CP(:,1), CP(:,2), porder, band);

            Eplot = interp2_matrix(x1d, y1d, Circle.Location(:,1), Circle.Location(:,2), porder, band);

    %         % G = make3DImplicitGaussian(y1d, x1d, z1d, sigma, spacing, band, 4, 1);
            GCart = make2DImplicitGaussian(x1d, y1d, sigma, spacing, band, numsigmas, LimitFarPoints);

            % GE = diag(G) + (G - diag(G))*Ecp;
            G = GCart*Ecp;

            save(fullfile(FileLocationCPGauss, FileNameEcp), 'Ecp', '-v7.3')
            save(fullfile(FileLocationCPGauss, FileNameEplot), 'Eplot', '-v7.3')
            save(fullfile(FileLocationCPGauss, FileNameGCart), 'GCart', '-v7.3')
            save(fullfile(FileLocationCPGauss, FileNameG), 'G', '-v7.3')
            
        else
            
            load(fullfile(FileLocationCPGauss, FileNameEcp))
            load(fullfile(FileLocationCPGauss, FileNameEplot))
            load(fullfile(FileLocationCPGauss, FileNameGCart))
            load(fullfile(FileLocationCPGauss, FileNameG))
            
        end
        
        
        % Gaussian Method
        % MinPoint = round(min(PointCloud.Location) - bandwidth - 2*spacing, 1);
        % MaxPoint = round(max(PointCloud.Location) + bandwidth + 2*spacing, 1);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Scale Parameter Estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ScaleParameter = findScaleParamter(sigma, alpha, NumSteps, 1, 2);
        
        % ScaleParameter = zeros(NumSteps,1);
        %
        % for i = 1 : NumSteps
        %
        %     ScaleParameter(i+1,1) = sqrt(i)*sigma;
        %
        % end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform Diffusion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Signal = zeros(length(CPSignal), NumSteps);
        % Signal(:,1) = CPSignal;
        Signal = CPSignal;
        
        % SignalAtVertex = zeros(Circle.LocationCount, NumSteps);
        % SignalAtVertex(:,1) = SignalOriginal;
        SignalAtVertex = SignalOriginal;
        
        WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumSteps-1));
        AbsErr = zeros(NumSteps,1);
        
        if ShowPlot
            figure(1)
        end
        for i = 1 : NumSteps - 1
            
            %     Signal(:,i+1) = GE * Signal(:,i);
            SignalNew = G * Signal;
            
            
            %     SignalAtVertex(:,i+1) = Eplot * Signal(:,i+1);
            SignalAtVertex = Eplot * SignalNew;
            
            Truth = ExactSignal( ScaleParameter(i+1,1), Theta);
            AbsErr(i+1,1) = norm(Truth - SignalAtVertex, inf);
            
            if ShowPlot 
                clf
                plot(Theta, SignalOriginal,'k')
                hold on
                plot(Theta, Truth, 'g-')
                plot(Theta, SignalAtVertex, 'r--')                
                legend('Original', 'Truth', 'Numerical')
            end
                
            Signal = SignalNew;
            
            waitbar(i/NumSteps, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumSteps-1));
        end
        
        waitbar(i/NumSteps, WaitBar, sprintf('Diffusion Complete'));
        close(WaitBar)
        if ShowPlot
            close(figure(1))
        end
        
        
        MCError(MCs, 1:2, MCsigma) = [NumSteps, AbsErr(NumSteps - 1)];
        MCErrorAll{MCs, MCsigma} = [(1:NumSteps)', AbsErr];
        
    end
end

close all



Points = zeros(size(MCErrorAll,1), 2, size(MCErrorAll,2));
for j = 1 : size(MCErrorAll,2)
    for i = 1 : size(MCErrorAll,1)
        Points(i,1:2,j) = [MCErrorAll{i,j}(end,1), MCErrorAll{i,j}(end,2)];
    end
end



colors = ['k','k','k','b','c','k','m','g','r'];
shapes = {'','','^','h','p','s','o','*','+'};
lin  = {':','-.','--','-','-','-','--','-.',':'};
msizes = [8, 8, 8, 8, 8, 17, 14, 12, 8];
%          [8, 8, 8, 8, 8, 10, 12, 14, 17];

figure('units','normalized','outerposition',[0 0 1 1])
for i = 1 : size(Points,3)
    loglog( Points(:,1,i), Points(:,2,i), strcat(colors(i), shapes{i}, lin{i}), 'linewidth', 3, 'markersize', msizes(i) )
    hold on
end

% % % 3rd Order Line
logx = [10,100];
logy = (10e-3).*logx.^(-3);
loglog(logx, logy,'k-','linewidth',3)

% xl=get(gca,'XLim');
% yl=get(gca,'YLim');
text1 = text(15,10*10^(-6),'$3^{rd}$ Order','Interpreter','latex');
set(text1,'Rotation',-20)
set(text1,'FontSize',50)



% % % 2nd Order Line
logx = [50,500];
logy = (10e1).*logx.^(-1);
loglog(logx, logy,'k-','linewidth',3)


text1 = text(80,10*10^(-0.4),'$1^{st}$ Order','Interpreter','latex');
set(text1,'Rotation',-8)
set(text1,'FontSize',50)


% % % Axis Labels
xlabel('N')
ylabel('$\| $error$ \|_{\infty}$','Interpreter','latex')
ax = gca;
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;


xlim([10*10^(-0.75) 10*10^(2.5)])

xticks([10e0 10e1 10e2 10e3])
yticks([10e-10, 10e-8, 10e-6, 10e-4, 10e-2, 10e0, 10e2, 10e4])
yticklabels({'10e-10','10e-8','10e-6','10e-4', '10e-2', '10e0', '10e2', '10e4'})

hleg = legend({'$m_\sigma = 1$','$m_\sigma = 2$','$m_\sigma = 3$','$m_\sigma = 4$','$m_\sigma = 5$','$m_\sigma = 6$','$m_\sigma = 7$','$m_\sigma = 8$','$m_\sigma = 9$'},'Interpreter','latex');
set(hleg, 'position', [0.80 0.58 0.07 0.2], 'FontSize', 35)













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

MCspacing = [0.25, 0.1, 0.05, 0.025, 0.01];

MCporder = [1,2,3,4,5,6,7,8];

MCError = zeros(length(MCspacing), 2, length(MCporder));
MCErrorAll = cell(length(MCspacing), length(MCporder));


for MCp = 1 : length(MCporder)
    
    for MCs = 1 : length(MCspacing)
        
        
        
        clearvars -except MCspacing MCs MCporder MCp MCError MCErrorAll ProjectRoot
        spacing = MCspacing(MCs)
        
       porder = MCporder(MCp)
        dim = 3;
%         spacing = 0.1;
        % sigma <= spacing
        sigma = spacing;
         numsigmas = 4
        LimitFarPoints = 1;
        
        ShowPlot = 0;
        
        if spacing > sigma
            bandwidth = 1.002*spacing*sqrt((dim-1)*((porder+1)/2)^2 + (((numsigmas*sigma)/spacing+(porder+1)/2)^2));
        else
            bandwidth = 1.002*numsigmas*sigma*sqrt((dim-1)*((porder+1)/2)^2 + (((numsigmas*sigma)/spacing+(porder+1)/2)^2));
        end
        
        NumberDivisions = 3;
        
        FileLocation = strcat(ProjectRoot, '/models/Sphere');
        FileName = strcat('Icosphere',num2str(NumberDivisions),'.off');
        
        
        tauImplicit = spacing / 8;
        MaxTauImplicit = 1/spacing;
        NumSteps = round(MaxTauImplicit); % 5000 ; %
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup File Name Directions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FileLocationCPGauss = strcat(ProjectRoot,'/models/Sphere/CPGauss/');
        FileNameCP = strcat('CP','_NumDiv',num2str(NumberDivisions),'_p',num2str(porder),'_s',num2str(spacing),'_NumSigma',num2str(numsigmas),'.mat');
        FileNameIJK = strcat('IJK','_NumDiv',num2str(NumberDivisions),'_p',num2str(porder),'_s',num2str(spacing),'_NumSigma',num2str(numsigmas),'.mat');
        
        FileNameEcp = strcat('Ecp','_NumDiv',num2str(NumberDivisions),'_p',num2str(porder),'_s',num2str(spacing),'_NumSigma',num2str(numsigmas),'.mat');
        FileNameEplot = strcat('Eplot','_NumDiv',num2str(NumberDivisions),'_p',num2str(porder),'_s',num2str(spacing),'_NumSigma',num2str(numsigmas),'.mat');
        FileNameGCart = strcat('GCart','_NumDiv',num2str(NumberDivisions),'_p',num2str(porder),'_s',num2str(spacing),'_NumSigma',num2str(numsigmas),'_LimitFar',num2str(LimitFarPoints),'.mat');
        FileNameG = strcat('G','_NumDiv',num2str(NumberDivisions),'_p',num2str(porder),'_s',num2str(spacing),'_NumSigma',num2str(numsigmas),'_LimitFar',num2str(LimitFarPoints),'.mat');
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make the Sphere and Signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ~exist(fullfile(FileLocation, FileName), 'file')
            
            [Location, Faces] = icosphere(NumberDivisions);
            [VerticesOut, FacesOut] = clearMeshDuplicates(Location, Faces );
            
            save_off(VerticesOut, FacesOut, fullfile(FileLocation, FileName))
            
            Sphere.Face = FacesOut;
            Sphere.Location = VerticesOut;
            
        else
            
            [Sphere.Location, Sphere.Face] = read_off( fullfile(FileLocation, FileName) );
            
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
        
        
        
        if ~exist(fullfile(FileLocation, FileName), 'file')
            save_off(Sphere.Location, Sphere.Face, fullfile(FileLocation, FileName))
        end
        
        
        [Theta, Phi, Radius] = cart2sph(Sphere.Location(:,1) ,Sphere.Location(:,2), Sphere.Location(:,3));
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Construct the Gaussian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        MinPoint = round(min(Sphere.Location) - bandwidth - ceil(numsigmas * sigma), 1);
        MaxPoint = round(max(Sphere.Location) + bandwidth + ceil(numsigmas * sigma), 1);
        
         
        
      %  if ~exist( fullfile(FileLocationCPGauss, FileNameCP), 'file') || ~exist( fullfile(FileLocationCPGauss, FileNameIJK), 'file')
            
            [IJK,DIST,CP,XYZ,CPFACE] = tri2cp(Sphere.Face, Sphere.Location, spacing, MinPoint, porder, 1);           
            
     %       save(fullfile(FileLocationCPGauss, FileNameCP), 'CP', '-v7.3')
    %        save(fullfile(FileLocationCPGauss, FileNameIJK), 'IJK', '-v7.3')
   %     else
  %          load(fullfile(FileLocationCPGauss, FileNameCP), 'CP')
 %           load(fullfile(FileLocationCPGauss, FileNameIJK), 'IJK')
%        end
        
        
        
        
        x1d = (MinPoint(1):spacing:MaxPoint(1))';
        y1d = (MinPoint(2):spacing:MaxPoint(2))';
        z1d = (MinPoint(3):spacing:MaxPoint(3))';
        
        BandSearchSize = [length(y1d), length(x1d), length(z1d)];
        
        Band = sub2ind(BandSearchSize, IJK(:,2), IJK(:,1), IJK(:,3));
        
        % FaceInterpolateWeights = interpBarycenterTriangle(Sphere, CP, CPFACE);
        
        % CPSignal = FaceInterpolateWeights * Sphere.Signal;
        

        
%        if ~exist( fullfile(FileLocationCPGauss, FileNameEcp), 'file')
            Ecp = interp3_matrix(x1d, y1d, z1d, CP(:,1), CP(:,2), CP(:,3), porder, Band);
%            save(fullfile(FileLocationCPGauss, FileNameEcp), 'Ecp', '-v7.3')
%        else
%            load(fullfile(FileLocationCPGauss, FileNameEcp))
%        end
        
%        if ~exist( fullfile(FileLocationCPGauss, FileNameEplot), 'file')
            Eplot = interp3_matrix(x1d, y1d, z1d, Sphere.Location(:,1), Sphere.Location(:,2), Sphere.Location(:,3), porder, Band);
%            save(fullfile(FileLocationCPGauss, FileNameEplot), 'Eplot', '-v7.3')
%        else
%            load(fullfile(FileLocationCPGauss, FileNameEplot))
%        end
        
        
%        clearvars -except MCspacing MCs MCsigma MCError MCErrorAll MCnumsigma ProjectRoot x1d y1d z1d spacing Band numsigmas LimitFarPoints sigma FileLocationCPGauss FileNameGCart FileNameEcp FileNameG

                
%        if ~exist( fullfile(FileLocationCPGauss, FileNameGCart), 'file')
            GCart  = make3DImplicitGaussian(x1d, y1d, z1d, sigma, spacing, Band, numsigmas, LimitFarPoints);
               
%            save(fullfile(FileLocationCPGauss, FileNameGCart), 'GCart', '-v7.3')
%        else
%            load(fullfile(FileLocationCPGauss, FileNameGCart))
%        end
        
   
        
%        if  ~exist( fullfile(FileLocationCPGauss, FileNameG), 'file')
            G = GCart*Ecp;
%            save(fullfile(FileLocationCPGauss, FileNameG), 'G', '-v7.3')
%        else
%            load(fullfile(FileLocationCPGauss, FileNameG))
%        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make Signals
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Define the signal on sphere
        ExactSignal = @(sigma, Phi) exp(-sigma^2)*sin(Phi);
        
        SignalOriginal = ExactSignal(0, Phi);
        
        Sphere.Signal = SignalOriginal;
        
        % Define the signal on CP of sphere
        [CPTheta, CPPhi, CPRadius] = cart2sph(CP(:,1), CP(:,2), CP(:,3));
        
        CPSignal = ExactSignal(0, CPPhi);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Scale Parameter Estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % ScaleParameter = findScaleParamter(sigma, alpha, NumSteps, 'natural', '2d');
        
        ScaleParameter = sqrt(0:NumSteps)'*sigma;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform Diffusion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Signal = zeros(length(CPSignal), NumSteps);
        % Signal(:,1) = CPSignal;
        Signal = CPSignal;
        
        Truth = zeros(Sphere.LocationCount, NumSteps);
        
%         SignalAtVertex = zeros(Sphere.LocationCount, NumSteps);
%         SignalAtVertex(:,1) = SignalOriginal;
        SignalAtVertex = SignalOriginal;
        
%        WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumSteps-1));
        AbsErr = zeros(NumSteps,1);
        
        if ShowPlot
            figure(1)
        end
        
        for i = 1 : NumSteps - 1
            
            SignalNew = G * Signal;
            
            SignalAtVertex = Eplot * SignalNew;
            
            Truth= ExactSignal(ScaleParameter(i+1), Phi);
            AbsErr(i+1,1) = norm(Truth - SignalAtVertex, inf);
            
            if ShowPlot
                clf
                plot(Phi, SignalOriginal,'k.')
                hold on
                plot(Phi, SignalAtVertex,'r.')
                plot(Phi, Truth, 'gd')
            end
            
            Signal = SignalNew;
            
 %           waitbar(i/NumSteps, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumSteps-1));
        end
        
%        waitbar(i/NumSteps, WaitBar, sprintf('Diffusion Complete'));
%        close(WaitBar)
        if ShowPlot
            close(figure(1))
        end
        
        MCError(MCs, 1:2, MCp) = [NumSteps, AbsErr(NumSteps - 1)];
        MCErrorAll{MCs, MCp} = [(1:NumSteps)', AbsErr];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
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

% min
% p = 3
% nsig = 4

% % % 0.5th Order Line
logx = [20,70];
logy = (10*10^(-3.2)).*logx.^(-0.5);
loglog(logx, logy,'k-','linewidth',3)

% xl=get(gca,'XLim');
% yl=get(gca,'YLim');
text1 = text(20,10*10^(-4.1),'$0.5^{th}$ Order','Interpreter','latex');
set(text1,'Rotation',-14)
set(text1,'FontSize',50)

% 
% 
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



xlim([10*10^(-0.75) 10*10^(1.5)])
% ylim([10*10^(-4) 10*10^(-2)])

xticks([10e0 10e1 10e2 10e3])
yticks([10e-10, 10e-8, 10e-6, 10e-5,10e-4,10e-3, 10e-2, 10e-1, 10e0, 10e2, 10e4])
yticklabels({'10e-10','10e-8','10e-6','10e-5','10e-4','10e-3', '10e-2', '10e-1', '10e0', '10e2', '10e4'})


hleg = legend({'$p = 1$','$p = 2$','$p = 3$','$p = 4$','$p = 5$','$p = 6$','$p = 7$','$p = 8$'},'Interpreter','latex');
set(hleg, 'position', [0.81 0.53 0.07 0.2], 'FontSize', 50)
% tit = title({'$m_\sigma = 4$'},'Interpreter','latex');
% set(tit, 'FontSize', 50)






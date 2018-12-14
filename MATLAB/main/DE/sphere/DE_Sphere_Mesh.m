% Andrew Rhodes
% ASEL
% March 2018


close all
clear
clc

global ProjectRoot;% Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BDF = 4;
ExactSignal = @(sigma, Phi) exp(-sigma^2)*sin(Phi);
FileLocation = strcat(ProjectRoot,'/models/sphere/');
FileLocationMesh = strcat(ProjectRoot,'/models/sphere/mesh/');
alpha = 1/2;
ShowPlot = 0;
options.dtype = 'geodesic';
options.htype = 'psp';
% options.dtype = 'euclidean';


MCdiv = [1, 2, 3, 4, 5, 6];
MCrho = [3, 4, 5, 6, 7];

MCError = zeros(length(MCdiv), 2, length(MCrho));
MCErrorAll = cell(length(MCdiv), length(MCrho));


for MCr = 1 : length(MCrho)
    
    for MCd = 1 : length(MCdiv)
        
        clearvars -except MCr MCd MCdiv MCrho MCErrorAll MCError ProjectRoot BDF ExactSignal FileLocation FileLocationMesh alpha ShowPlot options
        
        NumberDivisions = MCdiv(MCd)
        options.rho = MCrho(MCr)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup File Name Directions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        FileNameSphere = strcat('Icosphere/Icosphere',num2str(NumberDivisions),'.off');
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make the Sphere
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        if ~exist(fullfile(FileLocation, FileNameSphere), 'file')
            
            [Sphere.Location, Sphere.Face] = icosphere(NumberDivisions);
%             [VerticesOut, FacesOut] = clearMeshDuplicates(Location, Faces );
            
            save_off(Sphere.Location, phere.Face, fullfile(FileLocation, FileNameSphere))
            
        else
            
            [Sphere.Location, Sphere.Face] = read_off(fullfile(FileLocation, FileNameSphere));
            
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
        Sphere.FaceArea = findFaceArea(Sphere.Location, Sphere.Face);
        Sphere.Signal = zeros(Sphere.LocationCount,1);
        Sphere = findMeshResolution(Sphere, 'Model');
                
        
        options.hs = Sphere.Resolution/2;       
        
        % % % % % % % % % %
        tau = options.hs^2/4;
        MaxTau = 1 / Sphere.Resolution;
        NumSteps = round(MaxTau);
        % % % % % % % % % %
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make the Signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [Sphere.Theta, Sphere.Phi, Sphere.Radius] = cart2sph(Sphere.Location(:,1) ,Sphere.Location(:,2), Sphere.Location(:,3));
        
        
        SignalOriginal = ExactSignal(0, Sphere.Phi);
        Sphere.Signal = SignalOriginal;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Construct the Laplace Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FileNameItL = strcat(FileLocation, 'mesh/Icosphere_div',num2str(NumberDivisions),...
            '_ItL','_BDF',num2str(BDF),'_rho',num2str(options.rho),'_dtype_',options.dtype,'.mat');
        
        
        if ~exist(FileNameItL, 'file')
            ItL = makeMeshLaplaceBeltrami( fullfile(FileLocation, FileNameSphere), options, BDF, tau, alpha);
            save(FileNameItL, 'ItL', '-v7.3');
        else
            load(FileNameItL, 'ItL');
        end
                
 
        Signal = performBDFDiffusion(Sphere.Signal, NumSteps, ItL);
%         Signal = performBDFDiffusion_cpp(ItL, Sphere.Signal, NumSteps);
      

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Scale Parameter Estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        ScaleParameter = findScaleParameter(tau, alpha, NumSteps, 'Laplacian', 'Natural');

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform Diffusion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ShowPlot
            figure(1)
        end
        
        AbsErr = zeros(NumSteps,1);
        for i = 1 : NumSteps - 1
           
            Truth = ExactSignal(ScaleParameter(i+1), Sphere.Phi);
            AbsErr(i+1,1) = norm(Truth - Signal(:,i+1), inf);
            
            if ShowPlot
                clf
                plot(Sphere.Phi, SignalOriginal,'ko')
                hold on
                plot(Sphere.Phi, Truth, 'gd')
                plot(Sphere.Phi, Signal(:,i+1),'r.')
            end
            
        end
        
        if ShowPlot
            close(figure(1))
        end
        
        
        MCError(MCd, 1:2, MCr) = [NumSteps, AbsErr(NumSteps - 1)];
        MCErrorAll{MCd, MCr} = [(1:NumSteps)', AbsErr];
        
        
    end
    
end


save(strcat(FileLocation,'mesh/MCErrorAll.mat'), 'MCErrorAll', '-v7.3')


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












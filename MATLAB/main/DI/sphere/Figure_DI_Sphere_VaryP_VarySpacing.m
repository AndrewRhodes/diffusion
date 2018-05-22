% Andrew Rhodes
% ASEL
% February 2018

% Diffusion accuracy comparison for a signal on a sphere.


close all
clear
clc

global ProjectRoot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCspacing = [0.25, 0.1, 0.05, 0.025, 0.01, 0.005];%, 0.0025, 0.001];
MCporder = [1, 2, 3, 4, 5, 6, 7];

MCError = zeros(length(MCspacing), 2, length(MCporder));
MCErrorAll = cell(length(MCspacing), length(MCporder));


for MCp = 1 : length(MCporder)
    
    for MCs = 1 : length(MCspacing)
        
        clearvars -except MCspacing MCporder MCs MCp MCError MCErrorAll ProjectRoot
        spacing = MCspacing(MCs)
        porder = MCporder(MCp)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % User Defined Criteria
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        NumberDivisions = 5;
        alpha = 1;
        
%        porder = 5;
        dim = 3;
        Lorder = 4;
%         spacing = 0.05;
        Model = 'Icosphere';
        
        ShowPlot = 0;
        BDF = 2;
        tauFraction = 1/8;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup File Names
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FileLocationModel = strcat(ProjectRoot,'/models/sphere/');
        FileNameModelPly = strcat(Model,num2str(NumberDivisions),'.ply');
        FileNameModelOff = strcat(Model,num2str(NumberDivisions),'.off');
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make the Sphere and Signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        PointCloud = getIcosphere( fullfile(FileLocationModel, FileNameModelOff), NumberDivisions);
        
        
        
        % % % % % % % % % %
        tau = spacing * tauFraction;
        MaxTau = 1 / spacing;
        NumSteps = round(MaxTau);
        % % % % % % % % % %
        
        
        PointCloud.LocationCount = size(PointCloud.Location,1);
        PointCloud.FaceCount = size(PointCloud.Face, 1);
        PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
        
        
        [PointCloud.Theta, PointCloud.Phi, PointCloud.Radius] = cart2sph(PointCloud.Location(:,1) ,PointCloud.Location(:,2), PointCloud.Location(:,3));
        
        TrueSignalModel = @(sigma, Phi) exp(-sigma^2)*sin(Phi);
        
        
        SignalOriginal = TrueSignalModel(0, PointCloud.Phi);
        
        PointCloud.Signal = SignalOriginal;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup Laplace-Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        [ItL, Eplot, CP, CPFACE] = makeImplicitLaplaceBeltrami( PointCloud, porder, Lorder, spacing, dim, BDF, tau, alpha);
        
        FaceInterpolateWeights = interpBarycenterTriangle(PointCloud, CP, CPFACE);
        
        [CPTheta, CPPhi, CPRadius] = cart2sph(CP(:,1) ,CP(:,2), CP(:,3));
        
        CPSignal = TrueSignalModel(0, CPPhi);
        
        % CPSignal = FaceInterpolateWeights * PointCloud.Signal;
        
        
        
        ScaleParameter = findScaleParamter(tau, alpha, NumSteps, 'Laplacian', 'Natural');
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Construct the Laplace Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % BandSearchSize = [length(y1d), length(x1d), length(z1d)];
        % BandInit = sub2ind(BandSearchSize, IJK(:,1), IJK(:,2), IJK(:,3));
        
        % BandInit = find(abs(DIST) <= bandwidth);
        
        % FaceInterpolateWeights = interpBarycenterTriangle(Sphere, CP, CPFACE);
        
        % CPSignal = FaceInterpolateWeights * Sphere.Signal;
        
        % CPInit = CP(BandInit, :);
        %
        % XYZInit = XYZ(BandInit, :);
        
        
        % if ~exist(fullfile(FileLocationCP, FileNameL), 'file') || ~exist(fullfile(FileLocationCP, FileNameEcp), 'file') || ~exist(fullfile(FileLocationCP, FileNameEplot), 'file') || ~exist(fullfile(FileLocationCP, FileNameM), 'file')
        
        %     XYZInit = XYZ(BandInit,:);
        %     CPInit = CP(BandInit,:);
        %
        %     [L, Ecp, R, BandInner, BandOuter, BandInnerFull, BandOuterFull] = ...
        %         ops_and_bands3d(x1d, y1d, z1d, XYZInit(:,1), XYZInit(:,2), XYZInit(:,3), ...
        %         CPInit(:,1), CPInit(:,2), CPInit(:,3), BandInit, porder, Lorder);
        %
        %     CPIn = CPInit(BandInner,:);
        
        % Create L, E, M
        %     Ecp = interp3_matrix(x1d, y1d, z1d, CP(:,1), CP(:,2), CP(:,3), porder);
        %     [EcpRow, EcpCol, EcpVal] = find(Ecp);
        %     BandInner = unique(EcpCol);
        %
        %     L = laplacian_3d_matrix(x1d, y1d,z1d, Lorder, BandInner, BandInit);
        %     [LRow, LCol, LVal] = find(L);
        %     BandOuterTemp = unique(LCol);
        %     BandOuter = BandInit( BandOuterTemp );
        %
        %     Ecp = Ecp(BandOuterTemp, BandInner);
        %     L = L(:, BandOuterTemp);
        %
        %
        %
        %     InnerInOuter = zeros(size(BandInner));
        %     R = sparse([],[],[], length(BandInner), length(BandOuter), length(BandInner));
        %     for i = 1 : length(BandInner)
        %         I = find(BandOuter == BandInner(i));
        %         InnerInOuter(i) = I;
        %         R(i,I) = 1;
        %     end
        %
        %
        %     M = lapsharp_unordered(L, Ecp, R);
        %
        %     Eplot = interp3_matrix(x1d, y1d, z1d, Sphere.Location(:,1), Sphere.Location(:,2), Sphere.Location(:,3), porder, BandInnerFull);
        % %
        %     CPOut = CPInit(BandOuterTemp,:);
        %     CPIn = R*CPOut;
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Diffusion of Signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        CPSignalDiffused = performBDFDiffusion(CPSignal, NumSteps, ItL);
        
        Signal = performEplotProjection(CPSignalDiffused, Eplot);
        
        
        TrueSignal = makeTrueSignalSphere(TrueSignalModel, NumSteps, ScaleParameter, PointCloud.Phi);
        
        Error = findDiffusionError(TrueSignal, Signal, NumSteps, PointCloud.Phi, ShowPlot);
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        %         Ecp = interp3_matrix(x1d, y1d, z1d, CP(:,1), CP(:,2), CP(:,3), porder);
        %         [EcpRow, EcpCol, EcpVal] = find(Ecp);
        %         BandInner = unique(EcpCol);
        %
        %         L = laplacian_3d_matrix(x1d, y1d,z1d, Lorder, BandInner, BandInit);
        %         [LRow, LCol, LVal] = find(L);
        %         BandOuterTemp = unique(LCol);
        %         BandOuter = BandInit( BandOuterTemp );
        %
        %         CPOut = CPInit(BandOuterTemp,:);
        %
        % %         Reform the L, Ecp matrices
        %         Ecp = Ecp(BandOuterTemp, BandInner);
        %         L = L(:, BandOuterTemp);
        %
        %
        %
        %         [xp,yp,zp] = sphere(64);
        %         [th_plot, phi_plot, r] = cart2sph(xp(:),yp(:),zp(:));
        %
        %         SignalOriginal = ExactSignal(0, phi_plot);
        %         Sphere.Signal = SignalOriginal;
        %
        %
        %         Eplot = interp3_matrix(x1d, y1d, z1d, xp(:), yp(:), zp(:), porder, BandInner);
        %
        %         Sphere.Phi = phi_plot;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Restriction Operator
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %         InnerInOuter = zeros(size(BandInner));
        %         R = sparse([],[],[], length(BandInner), length(BandOuter), length(BandInner));
        %
        %         for i = 1 : length(BandInner)
        %            I = find(BandOuter == BandInner(i));
        %            InnerInOuter(i) = I;
        %            R(i,I) = 1;
        %         end
        %
        %         CPIn = R*CPOut;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define the Signal and Plot Matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %         [CPTheta, CPPhi, CPRadius] = cart2sph(CP(:,1), CP(:,2), CP(:,3));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Construct the Laplace-Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %         M = lapsharp_unordered(L, Ecp, R);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform Diffusion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         CPSignal = ExactSignal(0, CPPhi);
        %         CPSignal = CPSignal(BandInner,:);
        
        %         Signal = zeros(length(CPSignal), NumStepsImplicit);
        %         Signal(:,1) = CPSignal;
        
        %         SignalAtVertex = zeros(numel(xp), NumStepsImplicit);
        
        
        MCError(MCs, 1:2, MCp) = [NumSteps, Error(NumSteps - 1)];
        MCErrorAll{MCs, MCp} = [(1:NumSteps)', Error];
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
for i = [1,3,5,7]%1 : size(Points,3)
    loglog( Points(:,1,i), Points(:,2,i), strcat(colors(i), shapes{i}, lin{i}), 'linewidth', 3, 'markersize', msizes(i) )
    hold on
end




% % % 2nd Order Line
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






















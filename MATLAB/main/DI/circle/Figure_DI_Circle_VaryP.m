% Andrew Rhodes
% ASEL
% March 2018

% Diffusion accuracy comparison for a signal on a sphere.


close all
clear
clc

addpath('../../src/')
ProjectRoot = setupprojectpaths % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MCspacing = [1/3, 1/9, 1/22, 1/30, 1/90, 1/220, 1/300];
MCspacing = [0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001];
% MCporder = [1, 2, 3, 4, 5, 6, 7];
MCporder = [1, 3, 5, 7, 9];

MCError = zeros(length(MCspacing), 2, length(MCporder));
MCErrorAll = cell(length(MCspacing), length(MCporder));


for MCp = 1 : length(MCporder)
    
    for MCs = 1 : length(MCspacing)
        
        
        clearvars -except MCspacing MCporder MCs MCp MCError MCErrorAll ProjectRoot
        spacing = MCspacing(MCs)
        porder = MCporder(MCp)
        
        alpha = 1;
        dim = 2;
        Lorder = 4;
        NumCirclePoints = 500;
        ShowPlot = 0;
        
        bandwidth = 1.00001*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));
        
        tauImplicit = spacing / 8;
        MaxTauImplicit = 1/spacing;
        NumStepsImplicit = round(MaxTauImplicit); % 3000;%100*MaxTauImplicit; %ceil(MaxTauImplicit / tauImplicit);
        
        ExactSignal = @(sigma, theta) exp(-sigma^2/2)*cos(theta) + exp(-9*sigma^2/2)*sin(3*theta);
        % ExactSignal = @(sigma, theta) exp(-sigma^2/2)*cos(theta);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup File Name Directions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FileLocationCP = strcat(ProjectRoot, '/models/circle/CPLaplace/');
        FileNameCP = strcat('CP','_linspace',num2str(NumCirclePoints),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        FileNameDist = strcat('Dist','_linspace',num2str(NumCirclePoints),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        
        FileNameL = strcat('L','_linspace',num2str(NumCirclePoints),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        FileNameEplot = strcat('Eplot','_linspace',num2str(NumCirclePoints),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        FileNameEcp = strcat('Ecp','_linspace',num2str(NumCirclePoints),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        FileNameM = strcat('M','_linspace',num2str(NumCirclePoints),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        FileNameR = strcat('R','_linspace',num2str(NumCirclePoints),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make the Circle and Signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Theta = linspace(0, 2*pi, NumCirclePoints)';
        
        Radius = ones(size(Theta));
        [xp,yp] = pol2cart(Theta, Radius);
        Circle.Location(:,1) = xp(:);
        Circle.Location(:,2) = yp(:);
        Circle.LocationCount = length(Circle.Location);
        
        SignalOriginal = ExactSignal(0, Theta);
        Circle.Signal = SignalOriginal;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Build Grid and Inner/Outer Bands
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        MinPoint = min(Circle.Location) - bandwidth - spacing;
        MaxPoint = max(Circle.Location) + bandwidth + spacing;

        x1d = (MinPoint(1):spacing:MaxPoint(1))';
        y1d = (MinPoint(2):spacing:MaxPoint(2))';

        
        [GridX, GridY] = meshgrid(x1d, y1d);
        
%         if ~exist( fullfile(FileLocationCP,FileNameCP), 'file') || ~exist( fullfile(FileLocationCP, FileNameDist), 'file')
            [CP(:,1), CP(:,2), dist] = cpCircle(GridX(:), GridY(:));
%             save( fullfile(FileLocationCP, FileNameCP), 'CP', '-v7.3')
%             save( fullfile(FileLocationCP, FileNameDist), 'dist', '-v7.3')
%         else
            
%             load( fullfile(FileLocationCP, FileNameCP) )
%         end
        
        
        BandInit = find(abs(dist) <= bandwidth);

        CPInit = CP(BandInit,:);
        GridXBandInit = GridX(BandInit);
        GridYBandInit = GridY(BandInit);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Matric Construction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Construct full Ecp matric
        Ecp = interp2_matrix(x1d, y1d, CPInit(:,1), CPInit(:,2), porder);
        % Ecp = interp2_matrix(x1d, y1d, CP(:,1), CP(:,2), porder, BandInit);
        
        
        [EcpRow, EcpCol, EcpVal] = find(Ecp);
        BandInner = unique(EcpCol);
        
        % Construct full L matrix
        L = laplacian_2d_matrix(x1d, y1d, Lorder, BandInner, BandInit);
        
        [LRow, LCol, LVal] = find(L);
        BandOuterTemp = unique(LCol);
        BandOuter = BandInit( BandOuterTemp );
        
        
        CPOut = CPInit(BandOuterTemp,:);
        GridXOut = GridXBandInit(BandOuterTemp);
        GridYOut = GridYBandInit(BandOuterTemp);
        
        % Reform the L, Ecp matrices
        Ecp = Ecp(BandOuterTemp, BandInner);
        L = L(:, BandOuterTemp);
        
        
%         if ~exist( fullfile(FileLocationCP, FileNameEplot), 'file')
            Eplot = interp2_matrix(x1d, y1d, Circle.Location(:,1), Circle.Location(:,2), porder, BandInner);
%             save( fullfile(FileLocationCP, FileNameEplot), 'Eplot', '-v7.3')
%         else
%             load( fullfile(FileLocationCP, FileNameEplot) )
%         end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Restriction Operator
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         if ~exist( fullFile(FileLocationCP, FileNameR), 'file')
        
        InnerInOuter = zeros(size(BandInner));
        R = sparse([],[],[], length(BandInner), length(BandOuter), length(BandInner));
        
        for i = 1 : length(BandInner)
            I = find(BandOuter == BandInner(i));
            InnerInOuter(i) = I;
            R(i,I) = 1;
        end
        
        CPIn = R*CPOut;
        GridXIn = R*GridXOut;
        GridYIn = R*GridYOut;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define the Signal and Plot Matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        [CPTheta, CPr] = cart2pol(CPIn(:,1),CPIn(:,2));
        
        CPSignal = ExactSignal(0, CPTheta);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Construct the Laplace-Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ~exist( fullfile(FileLocationCP,FileNameM), 'file')
            M = lapsharp_unordered(L, Ecp, R);
            save( fullfile(FileLocationCP,FileNameM), 'M', '-v7.3')
        else
            load( fullfile(FileLocationCP,FileNameM) )
        end
        
        
        ItM = speye(size(M)) - alpha*tauImplicit * M;
        
        I23tM = speye(size(M)) - (2/3)*alpha*tauImplicit * M;
        
        I611M = speye(size(M)) - (6/11)*alpha*tauImplicit * M;
        
        I1225M = speye(size(M)) - (12/25)*alpha*tauImplicit * M;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Scale Parameter Estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ScaleParameter = findScaleParamter(tauImplicit, alpha, NumStepsImplicit, 'natural', '3d');
        
%         
%         ScaleParameter = zeros(NumStepsImplicit,1);
%         
%         for i = 1 : NumStepsImplicit - 1
%             
%             %     ScaleParameter(i+1,1) = sqrt(2*i*alpha*tauImplicit);
%             %     ScaleParameter(i+1,1) = alpha*i*tauImplicit;
%             ScaleParameter(i+1,1) = sqrt(2*i*alpha*tauImplicit);
%             
%         end
%         
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform Diffusion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Signal = zeros(length(CPSignal), NumStepsImplicit);
        Signal(:,1) = CPSignal;
        
        SignalAtVertex = zeros(Circle.LocationCount, NumStepsImplicit);
        SignalAtVertex(:,1) = SignalOriginal;
        
        WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsImplicit-1));
        AbsErr = zeros(NumStepsImplicit,1);
        
        NumIter = min( size(M,1), 100 );
        
        if ShowPlot
            figure(1)
        end
        
        for i = 1 : NumStepsImplicit - 1
            
                if i == 1
%                 Signal(:,i+1) = ItM \ Signal(:,i);
                    [Signal(:,i+1), flag] = bicg(ItM, Signal(:,i),  1e-10, NumIter);
                elseif i == 2
%                     Signal(:,i+1) = I23tM \ ( (4/3)*Signal(:,i) - (1/3)*Signal(:,i-1) );
                    [Signal(:,i+1), flag] = bicg(I23tM, (4/3)*Signal(:,i) - (1/3)*Signal(:,i-1), 1e-10, NumIter);
                elseif i == 3
% %                     Signal(:,i+1) = I611M \ ( (18/11)*Signal(:,i) - (9/11)*Signal(:,i-1) + (2/11)*Signal(:,i-2) );
                    [Signal(:,i+1), flag] = bicg(I611M, (18/11)*Signal(:,i) - (9/11)*Signal(:,i-1) + (2/11)*Signal(:,i-2), 1e-10, NumIter);
                else
%                     Signal(:,i+1) = I1225M \ ( (48/25)*Signal(:,i)-(36/25)*Signal(:,i-1) + (16/25)*Signal(:,i-2) - (3/25)*Signal(:,i-3) );
                    [Signal(:,i+1), flag] = bicg(I1225M, (48/25)*Signal(:,i) - (36/25)*Signal(:,i-1) + (16/25)*Signal(:,i-2) - (3/25)*Signal(:,i-3), 1e-10, NumIter);
                end
            
            if flag
                disp(flag)
            end
            
            % % Interpolate Signal on Circle
            SignalAtVertex(:,i+1) = Eplot * Signal(:,i+1);
            
            % % Calculate truth and find error
            Truth = ExactSignal(ScaleParameter(i+1), Theta);
            AbsErr(i+1,1) = norm(Truth - SignalAtVertex(:,i+1), inf);
            
            
            if ShowPlot
                clf
                plot(Theta, SignalOriginal,'b')
                hold on
                plot(Theta, SignalAtVertex(:,i+1),'k')
                plot(Theta, Truth,'r--')
%                 pause(0.1)
            end            
            
            
            
            waitbar(i/NumStepsImplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsImplicit-1));
            
        end
        
        waitbar(i/NumStepsImplicit, WaitBar, sprintf('Diffusion Complete'));
        close(WaitBar)
        
        if ShowPlot
            close(figure(1))
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot Errors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % SignalAtVertex(abs(SignalAtVertex)<eps) = 0;
        
        % AbsErr = zeros(NumStepsImplicit,1);
        % RelErr = zeros(NumStepsImplicit,1);
        % Truth = zeros(Circle.LocationCount,1);
        % Truth(:,1) = SignalOriginal;
        %
        % for i = 1 : NumStepsImplicit
        %
        %     Truth(:,i) = exp(-ScaleParameter(i,1)^2) * SignalOriginal;
        %
        % %     Truth(abs(Truth(:,i))<eps,i) = 0;
        %
        % %     Truth(:,i) = exp(-tauImplicit) * Truth(:,i-1);
        %
        % %     Truth(:,i) = exp(-ScaleParameter(i,1)^2) * cos(Theta) + exp(-9*ScaleParameter(i,1)^2) * sin(3*Theta);
        %
        %     AbsErr(i, 1) = norm( Truth(:,i) - SignalAtVertex(:,i), inf);
        %
        %
        % end
        
        % MCError(MC,1:2) = [NumStepsImplicit, AbsErr(NumStepsImplicit - 1)];
        % MCErrorAll{MC,1} = [(1:NumStepsImplicit)', AbsErr];
        
        MCError(MCs, 1:2, MCp) = [NumStepsImplicit, AbsErr(NumStepsImplicit - 1)];
        MCErrorAll{MCs, MCp} = [(1:NumStepsImplicit)', AbsErr];
        
    end
end


close all






Points = zeros(size(MCErrorAll,1), 2, size(MCErrorAll,2));
for j = 1 : size(MCErrorAll,2)
    for i = 1 : size(MCErrorAll,1)
        Points(i,1:2,j) = [MCErrorAll{i,j}(end,1), MCErrorAll{i,j}(end,2)];
    end
end

% FileLocationPoints = strcat(ProjectRoot,'/data/circle/CPLaplace/');
% FileNamePoints = strcat('Points','_NumPoints',num2str(NumCirclePoints),'_s',num2str(MCspacing(1)),'-',num2str(MCspacing(end)),'_p',num2str(MCporder(1)),'-',num2str(MCporder(end)),'_l',num2str(Lorder),'.mat');
% save(fullfile(FileLocationPoints, FileNamePoints), 'Points', '-v7.3')
% 
% load( fullfile(FileLocationPoints, FileNamePoints) )



colors = ['k','m','b','k','g','r'];
shapes = {'p','s','o','*','+'};
lin  = {'-','-','--','-.',':'};
msizes = [ 8, 17, 14, 12, 8];


figure('units','normalized','outerposition',[0 0 1 1])
for i = 1 : size(Points,3)
    loglog( Points(:,1,i), Points(:,2,i), strcat(colors(i), shapes{i}, lin{i}), 'linewidth', 3, 'markersize', msizes(i) )
%     pause
    hold on
end

% 2nd order line
logx = [10,100];
logy = (10e-3).*logx.^(-2);
loglog(logx, logy,'k-','linewidth',3)


text1 = text(10,10*10^(-6),'$2^{nd}$ Order','Interpreter','latex');
set(text1,'Rotation',-28)
set(text1,'FontSize',50)


% % % Axis Labels
xlabel('N')
ylabel('$\| $error$ \|_{\infty}$','Interpreter','latex')
ax = gca;
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;


xlim([10*10^(-0.75) 10*10^(3)])

xticks([10e0 10e1 10e2 10e3])
yticks([10e-10, 10e-8, 10e-6, 10e-4, 10e-2, 10e0, 10e2, 10e4])
yticklabels({'10e-10','10e-8','10e-6','10e-4', '10e-2', '10e0', '10e2', '10e4'})

hleg = legend({'$p = 1$','$p = 3$','$p = 5$','$p = 7$','$p = 9$'},'Interpreter','latex');
set(hleg, 'position', [0.78 0.56 0.07 0.2], 'FontSize', 70)


















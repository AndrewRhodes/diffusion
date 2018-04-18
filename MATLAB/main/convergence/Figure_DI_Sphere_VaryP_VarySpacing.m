% Andrew Rhodes
% ASEL
% February 2018

% Diffusion accuracy comparison for a signal on a sphere.


close all
clear
clc

addpath('../src/')
ProjectRoot = setupprojectpaths; % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCspacing = [0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025];%, 0.001];
MCporder = [1, 2, 3, 4, 5, 6, 7];

MCError = zeros(length(MCspacing), 2, length(MCporder));
MCErrorAll = cell(length(MCspacing), length(MCporder));


for MCp = 1 : length(MCporder)
    
    for MCs = 1 : length(MCspacing)
        
        clearvars -except MCspacing MCporder MCs MCp MCError MCErrorAll ProjectRoot
        spacing = MCspacing(MCs)
        porder = MCporder(MCp)
        
        NumberDivisions = 4;
        alpha = 1;
        
        MaxDegreeL = 50;
        
%         porder = 4; % order of interpolation
        dim = 3; % dimension
        Lorder = 2; % Cartesian Laplace order
%         spacing = 0.1; % spacing of embedding grid
        bandwidth = 1.002*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));
        
        tauImplicit = spacing / 8; % time step
        MaxTauImplicit = 1 / spacing;
        NumStepsImplicit = round(MaxTauImplicit); %ceil(MaxTauImplicit / tauImplicit);
        
        ShowPlot = 0;
        
        ExactSignal = @(sigma, Phi) exp(-sigma^2/2)*cos(bsxfun(@minus,Phi,pi/2));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup File Name Directions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FileLocation = strcat(ProjectRoot,'/models/Sphere/');
        FileName = strcat('Icosphere',num2str(NumberDivisions),'.off');
        
        FileLocationCP = strcat(ProjectRoot,'/models/Sphere/CPLaplace/');
        FileNameIJK = strcat('IJK','_Div',num2str(NumberDivisions),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        FileNameCP = strcat('CP','_Div',num2str(NumberDivisions),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        FileNameCPFACE = strcat('CPFACE','_Div',num2str(NumberDivisions),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        FileNameDIST = strcat('DIST','_Div',num2str(NumberDivisions),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        FileNameXYZ = strcat('XYZ','_Div',num2str(NumberDivisions),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        
        FileNameL = strcat('L','_Div',num2str(NumberDivisions),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        FileNameEplot = strcat('Eplot','_Div',num2str(NumberDivisions),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        FileNameEcp = strcat('Ecp','_Div',num2str(NumberDivisions),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        FileNameM = strcat('M','_Div',num2str(NumberDivisions),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        FileNameCPIn =strcat('CPIn','_Div',num2str(NumberDivisions),'_s',num2str(spacing),'_p',num2str(porder),'_l',num2str(Lorder),'.mat');
        
        
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
        
        
        %
        % [xp, yp, zp] = sphere(100);
        % Sphere.Location(:,1) = xp(:);
        % Sphere.Location(:,2) = yp(:);
        % Sphere.Location(:,3) = zp(:);
        Sphere.LocationCount = length(Sphere.Location);
        
        
        
        [Sphere.Theta, Sphere.Phi, Sphere.Radius] = cart2sph(Sphere.Location(:,1), Sphere.Location(:,2), Sphere.Location(:,3));
        
        % SphericalHarmonic = makeRealSphericalHarmonic( MaxDegreeL, Sphere.Theta, Sphere.Phi );
        
        
        % Define the signal
        % ExactSignal = @(sigma, SignalOriginal, MaxDegreeL) sum(cell2mat(cellfun(@times, num2cell( (exp(-(sigma^2/2).*(1:MaxDegreeL).*((1:MaxDegreeL)+1))')), SignalOriginal, 'UniformOutput', 0)),1);
        % ExactSignal = @(sigma, Theta) exp(-sigma^2/2)*cos(Theta);
        
        % SignalOriginal = ExactSignal(0, SphericalHarmonic, MaxDegreeL);
        SignalOriginal = ExactSignal(0, Sphere.Phi);
        % SignalOriginal = zeros(Sphere.LocationCount,1);
        % SignalPositions = find(Sphere.Theta == 0 & Sphere.Phi == pi/2);
        % SignalOriginal(SignalPositions) = ones(size(SignalPositions));
        % SignalOriginal(1) = 1;
        Sphere.Signal = SignalOriginal;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Construct the Laplace Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        MinPoint = min(Sphere.Location) - bandwidth - spacing;
        MaxPoint = max(Sphere.Location) + bandwidth + spacing;
        
        % MinPoint = [-3, -3, -3];
        % MaxPoint = [3, 3, 3];
        
        x1d = (MinPoint(1):spacing:MaxPoint(1))';
        y1d = (MinPoint(2):spacing:MaxPoint(2))';
        z1d = (MinPoint(3):spacing:MaxPoint(3))';
        
        
%         [IJK,DIST,CP,XYZ,CPFACE] = tri2cp(Sphere.Face, Sphere.Location, spacing, MinPoint, porder, Lorder/2);
        
        
%         BandSearchSize = [length(x1d), length(y1d), length(z1d)];
%         
%         Band = sub2ind(BandSearchSize, IJK(:,1), IJK(:,2), IJK(:,3));
        
%         BandInit = find(abs(DIST) <= 1.5*bandwidth);
        
%         FaceInterpolateWeights = interpBarycenterTriangle(Sphere, CP, CPFACE);
        
%         CPSignal = FaceInterpolateWeights * Sphere.Signal;
        % CPSignal = zeros(length(CP),1);
        % CPSignalLocations = find(CP(:,1) == -1 & CP(:,2) == 0 & CP(:,3) == 0);
        % CPSignal(CPSignalLocations, 1) = ones(size(CPSignalLocations));
        % CPSignal(CPSignalLocations(1), 1) = 1;
        
        [GridX, GridY, GridZ] = meshgrid(x1d, y1d, z1d);
        [CP(:,1), CP(:,2), CP(:,3), dist] = cpSphere(GridX(:), GridY(:), GridZ(:));
        
        %
        %
        BandInit = find(abs(dist) <= bandwidth);
        %
        CPInit = CP(BandInit, :);
        
        [th, phi, r] = cart2sph(GridX, GridY, GridZ);

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Matric Construction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         L = laplacian_3d_matrix(x1d, y1d, z1d, Lorder, Band);
        
%         Eplot = interp3_matrix(x1d, y1d, z1d, Sphere.Location(:,1), Sphere.Location(:,2), Sphere.Location(:,3), porder, Band);
%         % Eplot = interpLagrange3D(BandSearchSize, MinPoint, PointCloud.Location, porder, Band, spacing);
        
%         Ecp = interp3_matrix(x1d, y1d, z1d, CP(:,1), CP(:,2), CP(:,3), porder, Band);
%         % Ecp = interpLagrange3D(BandSearchSize, MinPoint, CP, porder, Band, spacing);
        
%         M = lapsharp(L, Ecp);
        
%         M = M - diag(diag(M));
%         M = M - diag(sum(M,2));
        
        
        
        Ecp = interp3_matrix(x1d, y1d, z1d, CP(:,1), CP(:,2), CP(:,3), porder);
        [EcpRow, EcpCol, EcpVal] = find(Ecp);
        BandInner = unique(EcpCol);
        
        L = laplacian_3d_matrix(x1d, y1d,z1d, Lorder, BandInner, BandInit);
        [LRow, LCol, LVal] = find(L);
        BandOuterTemp = unique(LCol);
        BandOuter = BandInit( BandOuterTemp );
        
        CPOut = CPInit(BandOuterTemp,:);
        
%         Reform the L, Ecp matrices
        Ecp = Ecp(BandOuterTemp, BandInner);
        L = L(:, BandOuterTemp);
        
        
        
        [xp,yp,zp] = sphere(64);
        [th_plot, phi_plot, r] = cart2sph(xp(:),yp(:),zp(:));
        
        SignalOriginal = ExactSignal(0, phi_plot);
        Sphere.Signal = SignalOriginal;
        
        
        Eplot = interp3_matrix(x1d, y1d, z1d, xp(:), yp(:), zp(:), porder, BandInner);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Restriction Operator
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        InnerInOuter = zeros(size(BandInner));
        R = sparse([],[],[], length(BandInner), length(BandOuter), length(BandInner));
        
        for i = 1 : length(BandInner)
           I = find(BandOuter == BandInner(i));
           InnerInOuter(i) = I;
           R(i,I) = 1;
        end
        
        CPIn = R*CPOut;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define the Signal and Plot Matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [CPTheta, CPPhi, CPRadius] = cart2sph(CP(:,1), CP(:,2), CP(:,3));
        
        
        % [CPTheta, CPPhi, CPRadius] = cart2sph(CPIn(:,1), CPIn(:,2), CPIn(:,3));
        
        % SphericalHarmonicCP = makeRealSphericalHarmonic( MaxDegreeL, CPTheta, CPPhi );
        
        % CPSignal = ExactSignal(0, SphericalHarmonicCP, MaxDegreeL);
        % CPSignal = ExactSignal(0, Sphere.Theta);
        % CPSignal = zeros(length(CPTheta),1);
        % CpSignalPositions = find(CPPhi == 0);
        % CPSignal(CpSignalPositions) = ones(length(CpSignalPositions),1);
        % CPSignal(1) = 1;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Construct the Laplace-Beltrami
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        M = lapsharp_unordered(L, Ecp, R);
        
        
        ItM = speye(size(M)) - alpha*tauImplicit * M;
        
        I23tM = speye(size(M)) - (2/3)*alpha*tauImplicit * M;
        
        I611tM = speye(size(M)) - (6/11)*alpha*tauImplicit * M;
        
        I1225tM = speye(size(M)) - (12/25)*alpha*tauImplicit * M;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Scale Parameter Estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ScaleParameter = findScaleParamter(tauImplicit, 2*alpha, NumStepsImplicit, 'natural', '3d');
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform Diffusion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        CPSignal = ExactSignal(0, CPPhi);
        CPSignal = CPSignal(BandInner,:);
        
        Signal = zeros(length(CPSignal), NumStepsImplicit);
        Signal(:,1) = CPSignal;
        
%         SignalAtVertex = zeros(Sphere.LocationCount, NumStepsImplicit);
        SignalAtVertex = zeros(numel(xp), NumStepsImplicit);
        
        SignalAtVertex(:,1) = SignalOriginal;
        
        AbsErr = zeros(NumStepsImplicit,1);
        
        
        WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsImplicit-1));
        if ShowPlot
            figure(1)
        end
        
        for i = 1 : NumStepsImplicit - 1
            
                if i == 1
                    [Signal(:,i+1), flag] = bicg(ItM, Signal(:,i), 1e-10, 100);
                else%if i == 2
                    [Signal(:,i+1), flag] = bicg(I23tM, (4/3)*Signal(:,i) - (1/3)*Signal(:,i-1), 1e-10, 100);
            %     elseif i == 3
            %         [Signal(:,i+1), flag] = bicg(I611tM, (18/11)*Signal(:,i) - (9/11)*Signal(:,i-1) + (2/11)*Signal(:,i-2), 1e-10, 100);
            %     else
            %         [Signal(:,i+1), flag] = bicg(I1225tM, (48/25)*Signal(:,i)-(36/25)*Signal(:,i-1) + (16/25)*Signal(:,i-2) - (3/25)*Signal(:,i-3), 1e-10, 100);
                end
            
            
            
            % Interpolate back to explicit surface
            SignalAtVertex(:,i+1) = Eplot * Signal(:,i+1);
            
            if flag
                disp(flag)
            end
            
            % Calculate Truth and Error
            %     Truth = ExactSignal(ScaleParameter(i), SphericalHarmonic, MaxDegreeL);
            Truth = ExactSignal(ScaleParameter(i+1), phi_plot);
            
            AbsErr(i+1,1) = norm(Truth - SignalAtVertex(:,i+1), inf);
            
            
            if ShowPlot
                clf
                plot(Sphere.Phi, SignalOriginal,'ko')
                hold on
                plot(Sphere.Phi, Truth, 'gd')
                plot(Sphere.Phi, SignalAtVertex(:,i),'r.')
                legend('Original', 'Truth at i', 'Diffused at i')
                
            end
            
            waitbar(i/NumStepsImplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsImplicit-1));
            
        end
        
        waitbar(i/NumStepsImplicit, WaitBar, sprintf('Diffusion Complete'));
        close(WaitBar)
        
        if ShowPlot
            close(figure(1))
        end
        
        
        MCError(MCs, 1:2, MCp) = [NumStepsImplicit, AbsErr(NumStepsImplicit - 1)];
        MCErrorAll{MCs, MCp} = [(1:NumStepsImplicit)', AbsErr];
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure
loglog((1:NumStepsImplicit)', AbsErr)
xlabel('Iteration Number')
ylabel('Relative Error')

figure
plot(1:NumStepsImplicit, AbsErr)
xlabel('Iteration Number')
ylabel('Relative Error')



















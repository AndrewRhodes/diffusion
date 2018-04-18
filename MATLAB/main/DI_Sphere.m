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


NumberDivisions = 3;
alpha = 1;

porder = 5; 
dim = 3;
Lorder = 2;
spacing = 0.1;
bandwidth = 1.002*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));

tauImplicit = spacing / 8;
MaxTauImplicit = 1/spacing;
NumStepsImplicit = round(MaxTauImplicit); %ceil(MaxTauImplicit / tauImplicit);

ShowPlot = 1;

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
    
    [Sphere.Location, Sphere.Face] = icosphere(NumberDivisions);
    
    save_off(Sphere.Location, Sphere.Face, fullfile(FileLocation, FileName))
    
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

MinPoint = min(Sphere.Location) - bandwidth - spacing;
MaxPoint = max(Sphere.Location) + bandwidth + spacing;

x1d = (MinPoint(1):spacing:MaxPoint(1))';
y1d = (MinPoint(2):spacing:MaxPoint(2))';
z1d = (MinPoint(3):spacing:MaxPoint(3))';


% [GridX, GridY, GridZ] = meshgrid(x1d, y1d, z1d);
% [CP(:,1), CP(:,2), CP(:,3), DIST] = cpSphere(GridX(:), GridY(:), GridZ(:));


[Theta, Phi, Radius] = cart2sph(Sphere.Location(:,1) ,Sphere.Location(:,2), Sphere.Location(:,3));


SignalOriginal = ExactSignal(0, Phi);
% SignalOriginal = zeros(Sphere.LocationCount,1);
% SignalOriginal(1) = 1;
Sphere.Signal = SignalOriginal;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the Laplace Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% if ~exist(fullfile(FileLocationCP,FileNameIJK), 'file') || ~exist(fullfile(FileLocationCP,FileNameCP), 'file') || ~exist(fullfile(FileLocationCP,FileNameCPFACE), 'file') || ~exist(fullfile(FileLocationCP,FileNameDIST), 'file')
%     
%     [IJK,DIST,CP,XYZ,CPFACE] = tri2cp(Sphere.Face, Sphere.Location, spacing, MinPoint, porder, Lorder/2);
%     
%     save( fullfile(FileLocationCP, FileNameIJK), 'IJK', '-v7.3')
%     save( fullfile(FileLocationCP, FileNameCP), 'CP', '-v7.3')
%     save( fullfile(FileLocationCP, FileNameCPFACE), 'CPFACE', '-v7.3')
%     save( fullfile(FileLocationCP, FileNameDIST), 'DIST', '-v7.3')
%     save( fullfile(FileLocationCP, FileNameXYZ), 'XYZ', '-v7.3')
%     
% else
%     
%     load( fullfile(FileLocationCP, FileNameIJK) )
%     load( fullfile(FileLocationCP, FileNameCP) )
%     load( fullfile(FileLocationCP, FileNameCPFACE) )
%     load( fullfile(FileLocationCP, FileNameDIST) )
%     load( fullfile(FileLocationCP, FileNameXYZ) )
% end




% BandSearchSize = [length(y1d), length(x1d), length(z1d)];
% BandInit = sub2ind(BandSearchSize, IJK(:,1), IJK(:,2), IJK(:,3));

BandInit = find(abs(DIST) <= bandwidth);

% FaceInterpolateWeights = interpBarycenterTriangle(Sphere, CP, CPFACE);

% CPSignal = FaceInterpolateWeights * Sphere.Signal;

% CPInit = CP(BandInit, :);
% 
% XYZInit = XYZ(BandInit, :);


% if ~exist(fullfile(FileLocationCP, FileNameL), 'file') || ~exist(fullfile(FileLocationCP, FileNameEcp), 'file') || ~exist(fullfile(FileLocationCP, FileNameEplot), 'file') || ~exist(fullfile(FileLocationCP, FileNameM), 'file')
    
    XYZInit = XYZ(BandInit,:);
    CPInit = CP(BandInit,:);
    
    [L, Ecp, R, BandInner, BandOuter, BandInnerFull, BandOuterFull] = ...
        ops_and_bands3d(x1d, y1d, z1d, XYZInit(:,1), XYZInit(:,2), XYZInit(:,3), ...
        CPInit(:,1), CPInit(:,2), CPInit(:,3), BandInit, porder, Lorder);
    
    CPIn = CPInit(BandInner,:);  
    
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
    
    
%     save( fullfile(FileLocationCP, FileNameL), 'L', '-v7.3')
%     save( fullfile(FileLocationCP, FileNameEplot), 'Eplot', '-v7.3')
%     save( fullfile(FileLocationCP, FileNameEcp), 'Ecp', '-v7.3')
%     save( fullfile(FileLocationCP, FileNameM), 'M', '-v7.3')
%     save( fullfile(FileLocationCP, FileNameCPIn), 'CPin', '-v7.3')
%     
% else
%     
%     load( fullfile(FileLocationCP, FileNameL) )
%     load( fullfile(FileLocationCP, FileNameEplot) )
%     load( fullfile(FileLocationCP, FileNameEcp) )
%     load( fullfile(FileLocationCP, FileNameM) )
%     load( fullfile(FileLocationCP, FileNameCPIn) )
%     
% end

ItM = speye(size(M)) - alpha*tauImplicit * M;

I23tM = speye(size(M)) - (2/3)*alpha*tauImplicit * M;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the Signal and Plot Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[CPTheta, CPPhi, CPRadius] = cart2sph(CPIn(:,1) ,CPIn(:,2), CPIn(:,3));

% CPSignal = ExactSignal(0, CPTheta);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameter Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ScaleParameter = zeros(NumStepsImplicit,1);

for i = 1 : NumStepsImplicit - 1
   
%     ScaleParameter(i+1,1) = sqrt(2*alpha*i*tauImplicit);
    ScaleParameter(i+1,1) = sqrt(2*alpha*i*tauImplicit);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Signal = zeros(length(CPSignal), NumStepsImplicit);
Signal(:,1) = CPSignal;

SignalAtVertex = zeros(Sphere.LocationCount, NumStepsImplicit);
SignalAtVertex(:,1) = SignalOriginal;

AbsErr = zeros(NumStepsImplicit,1);


WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsImplicit-1));
figure(1)
for i = 1 : NumStepsImplicit - 1
    
%     if i == 1
    Signal(:,i+1) = ItM \ Signal(:,i);
%         [Signal(:,i+1), flag] = gmres(ItM, Signal(:,i), [], 1e-10, 100);
%     else
%         [Signal(:,i+1), flag] = gmres(I23tM, (4/3)*Signal(:,i) - (1/3)*Signal(:,i-1), [], 1e-10, 100);
%     end
    
    SignalAtVertex(:,i+1) = Eplot * Signal(:,i+1);
    
    if flag
        disp(flag)
    end
    
    Truth = ExactSignal(ScaleParameter(i+1), Phi);
    AbsErr(i,1) = norm(Truth - SignalAtVertex(:,i), inf);
    
    if ShowPlot
        clf
        plot(Phi, SignalOriginal,'b.')
        hold on
        plot(Phi, SignalAtVertex(:,i),'ko')  
        plot(Phi, Truth,'rs')
    end
	
   
    waitbar(i/NumStepsImplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsImplicit-1));
    
end

waitbar(i/NumStepsImplicit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)
close(figure(1))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



AbsErr = zeros(NumStepsImplicit,1);
RelErr = zeros(NumStepsImplicit,1);
Truth = zeros(Sphere.LocationCount,1);
Truth(:,1) = SignalOriginal;

for i = 2 : NumStepsImplicit 
    
    Truth(:,i) = exp(-2*ScaleParameter(i,1)) * SignalOriginal;

%     Truth(:,i) = exp(-tauImplicit) * Truth(:,i-1);

%     Truth(:,i) = exp(-(1.5^2)*ScaleParameter(i,1)) * cos(1.5*Phi) + exp(-(3^2)*ScaleParameter(i,1)) * sin(3*Theta);

    AbsErr(i, 1) = norm( Truth(:,i) - SignalAtVertex(:,i), inf);
    
% %     RelErr(i, 1) = norm( AbsErr(i, 1), inf) / norm( Truth(:,i), inf);
% %     RelErr(i, 1) = norm( (Truth(:,i) - Signal(:,i+1)) ./ Truth(:,i), inf);
    
end


% figure
loglog(1:NumStepsImplicit, AbsErr)
xlabel('Iteration Number')
ylabel('Relative Error')

hold on



















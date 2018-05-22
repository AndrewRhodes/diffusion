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


NumberDivisions = 3;
alpha = 1;

porder = 5; 
dim = 3;
Lorder = 2;
spacing = 0.05;
Model = 'Icosphere';

ShowPlot = 1;
BDF = 1;
tauFraction = 1/2;



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















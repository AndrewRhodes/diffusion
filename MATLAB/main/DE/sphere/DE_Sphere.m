% Andrew Rhodes
% ASEL
% March 2018


close all
clear
clc

global ProjectRoot; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

alpha = 1;
NumberDivisions = 4; % For building Icosphere

options.rho = 6;
options.dtype = 'geodesic';
% options.dtype = 'euclidean';

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
% Load the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PointCloud = getIcosphere( fullfile(FileLocationModel, FileNameModelOff), NumberDivisions);


PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
PointCloud = findMeshResolution(PointCloud, 'Model');


% % % % % % % % % % 
tau = PointCloud.Resolution * tauFraction;
MaxTau = 1 / PointCloud.Resolution;
NumSteps = round(MaxTau);
% % % % % % % % % % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[PointCloud.Theta, PointCloud.Phi, PointCloud.Radius] = cart2sph(PointCloud.Location(:,1) ,PointCloud.Location(:,2), PointCloud.Location(:,3));


TrueSignalModel = @(sigma, Phi) exp(-sigma^2 )*sin(Phi);


SignalOriginal = TrueSignalModel(0, PointCloud.Phi);

PointCloud.Signal = SignalOriginal;


% 
% DistanceFromImpulse = acos( Sphere.Location(1,:) * Sphere.Location')';
% EuclidDistanceFromImpulse = sqrt(sum(bsxfun(@minus, Sphere.Location, Sphere.Location(1,:)).^2,2));
% DistanceFromImpulse = asin( (EuclidDistanceFromImpulse/2) .* sqrt(4-EuclidDistanceFromImpulse.^2));
% 
% 
% ExactSignal = @(sigma, OriginalSignal)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ItL = makeExplicitLaplaceBeltrami( fullfile( FileLocationModel, FileNameModelOff ), options, BDF, tau, alpha);


ScaleParameter = findScaleParamter(tau, alpha, NumSteps, 'Laplacian', 'Natural');








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusion of Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Signal = performBDFDiffusion(PointCloud.Signal, NumSteps, ItL);

TrueSignal = makeTrueSignalSphere(TrueSignalModel, NumSteps, ScaleParameter, PointCloud.Phi);

Error = findDiffusionError(TrueSignal, Signal, NumSteps, PointCloud.Phi, ShowPlot);






figure
loglog(1:NumSteps, Error)



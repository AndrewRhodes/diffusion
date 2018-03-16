% Andrew Rhodes
% ASEL
% March 2018


close all
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional Paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd ~/Desktop/Ashish/'CS3 Code'/
addpath(genpath('~/Documents/Software/MeshLP/'))
addpath('~/AFOSR/Ashish/CS3 Code/')
addpath(genpath('~/GitProjects/pose/MATLAB_PointCloudDescriptors/OURCVFH/'))
addpath(genpath('../'))
addpath('../src/')
addpath('../data/')
addpath(genpath('../models/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCdiv = [3, 4, 5, 6, 7, 8]
ErrorAll = zeros
MCErrorAll{MCs, MCp} = [(1:NumStepsImplicit)', AbsErr];

for MC = 1 : length(MCdiv)
    clearvars -except MC MCdiv 

    NumberDivisions = MCdiv(MC);

    
    
FileLocation = '../models/Sphere/';
FileName = strcat('Icosphere',num2str(NumberDivisions),'.off');

alpha = 1;


options.rho = 5;
options.dtype = 'geodesic';
% options.dtype = 'euclidean';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Sphere and Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sphere.Location,Sphere.Face] = icosphere(NumberDivisions);

Sphere.FaceCount = size(Sphere.Face, 1);
Sphere.LocationCount = size(Sphere.Location,1);
Sphere.FaceArea = findFaceAreas(Sphere.Location,Sphere.Face);

Sphere.FaceArea = findFaceArea(Sphere.Location, Sphere.Face);
Sphere.Signal = zeros(Sphere.LocationCount,1);
Sphere = findMeshResolution(Sphere, 'Model');


% % % % % % % % % % 
tauExplicit = Sphere.Resoltution / 8;
MaxTauExplicit = 1 / Sphere.Resoltution;
NumStepsExplicit = round(MaxTauExplicit);
% % % % % % % % % % 


if ~exist(fullfile(FileLocation, FileName), 'file')
    save_off(Sphere.Location, Sphere.Face, fullfile(FileLocation, FileName))
end


[Theta, Phi, Radius] = cart2sph(Sphere.Location(:,1) ,Sphere.Location(:,2), Sphere.Location(:,3));

% Define the signal
% SignalOriginal = cos(1.5*Phi + pi/2) + sin(2.5*Theta - pi/2);
SignalOriginal = cos(Theta);% + sin(2.5*Theta);


Sphere.Signal = SignalOriginal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the Laplace Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationMeshLP = '../models/Sphere/meshLP/';
FileNameLapMat = strcat('LapMatMeshWeights','_eSS',num2str(Sphere.Resolution),'_',num2str(Sphere.LocationCount),'.mat');
FileNameArea = strcat('Area','_eSS',num2str(Sphere.Resolution),'_',num2str(Sphere.LocationCount),'.mat');
FileNamehEdge = strcat('hEdge2','_eSS',num2str(Sphere.Resoltution),'_',num2str(Sphere.LocationCount),'.mat');

if ~exist(fullfile(FileLocationMeshLP, FileNameLapMat), 'file')
    
    [LapMatMeshWeights, Area, hEdge2] = symmshlp_matrix(fullfile(FileLocation, FileName), options);
    
    save(fullfile(FileLocationMeshLP, FileNameLapMat), 'LapMatMeshWeights')
    save(fullfile(FileLocationMeshLP, FileNameArea), 'Area')
    save(fullfile(FileLocationMeshLP, FileNamehEdge), 'hEdge2')
else
    
    load(fullfile(FileLocationMeshLP, FileNameLapMat))
    load( fullfile(FileLocationMeshLP, FileNameArea) )
    load( fullfile(FileLocationMeshLP, FileNamehEdge) )
    
end

hEdge = (hEdge2/2);

Alength = length(A);

A1 = sparse(1:Alength, 1:Alength, 1./Area);

LBM = A1 * LapMatMeshWeights;

ItL = speye(Alength, Alength) - alpha*tauExplicit * LBM;

I23tL = speye(Alength, Alength) - (2/3)*alpha*tauExplicit * LBM;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameter Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ScaleParameter = findScaleParamter(tauExplicit, alpha, NumStepsExplicit, 1, 3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SignalExplicit = zeros(PointCloud.LocationCount, NumStepsExplicit);
SignalExplicit(:,1) = PointCloud.Signal;

AbsErr = zeros(NumStepsExplicit,1);


WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsExplicit-1));

for i = 1 : NumStepsExplicit - 1
    
    if i ==1
    %     [SignalExplicit(:,i+1), flag] = bicg(ItL, SignalExplicit(:,i), [], 1e-10, 60);
    [SignalExplicit(:,i+1), flag] = gmres(ItL, SignalExplicit(:,i), [], 1e-10, 100);
    else
        [SignalExplicit(:,i+1), flag] = gmres(I23tL, (4/3)*SignalExplicit(:,i) - (1/3)*SignalExplicit(:,i-1), [], 1e-10, 100);
    end
    
    if flag
        flag
    end
    
    Truth = exp(-ScaleParameter(i+1)^2/2) .* SignalOriginal;
    
	AbsErr(i+1,1) = norm(Truth - SignalExplicit(:,i+1), inf);
    
    
    waitbar(i/NumStepsExplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsExplicit-1));
    
end

waitbar(i/NumStepsExplicit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)



MCError(MCs, 1:2, MCp) = [NumStepsImplicit, AbsErr(NumStepsImplicit - 1)];
MCErrorAll{MCs, MCp} = [(1:NumStepsImplicit)', AbsErr];














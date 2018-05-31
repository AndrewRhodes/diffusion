% Andrew Rhodes
% ASEL
% March 2018

% Diffusion accuracy comparison for a signal on a sphere.


close all
clear
clc

addpath('../src/')
ProjectRoot = setupprojectpaths; % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
porder = 5; 
dim = 3;
spacing = 0.1;
% sigma <= spacing
sigma = spacing;
numsigmas = 7;
LimitFarPoints = 1;

ShowPlot = 1;

if spacing > sigma
    bandwidth = 1.002*spacing*sqrt((dim-1)*((porder+1)/2)^2 + (((numsigmas*sigma)/spacing+(porder+1)/2)^2));
else
    bandwidth = 1.002*numsigmas*sigma*sqrt((dim-1)*((porder+1)/2)^2 + (((numsigmas*sigma)/spacing+(porder+1)/2)^2));
end

NumberDivisions = 4;

FileLocation = '../models/Sphere/';
FileName = strcat('Icosphere',num2str(NumberDivisions),'.off');


MaxTauImplicit = 1/spacing;
NumSteps = round(MaxTauImplicit); % 5000 ; %

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

% Define the signal
ExactSignal = @(sigma, Phi) exp(-sigma^2/2)*cos(bsxfun(@minus,Phi,pi/2));

SignalOriginal = ExactSignal(0, Phi);

Sphere.Signal = SignalOriginal;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MinPoint = round(min(Sphere.Location) - bandwidth - ceil(numsigmas * sigma), 1);
MaxPoint = round(max(Sphere.Location) + bandwidth + ceil(numsigmas * sigma), 1);

% MinPoint = round(min(Circle.Location(:,1)) - bandwidth, 1);
% MaxPoint = round(max(Circle.Location(:,1)) + bandwidth, 1);


[IJK,DIST,CP,XYZ,CPFACE] = tri2cp(Sphere.Face, Sphere.Location, spacing, MinPoint, porder, 1);

x1d = (MinPoint(1):spacing:MaxPoint(1))';
y1d = (MinPoint(2):spacing:MaxPoint(2))';
z1d = (MinPoint(3):spacing:MaxPoint(3))';

BandSearchSize = [length(y1d), length(x1d), length(z1d)];

Band = sub2ind(BandSearchSize, IJK(:,2), IJK(:,1), IJK(:,3));

% FaceInterpolateWeights = interpBarycenterTriangle(Sphere, CP, CPFACE);

% CPSignal = FaceInterpolateWeights * Sphere.Signal;

[CPTheta, CPPhi, CPRadius] = cart2sph(CP(:,1), CP(:,2), CP(:,3));

CPSignal = ExactSignal(0, CPPhi);

% L = laplacian_3d_matrix(x1d,y1d,z1d, Lorder, Band);

Ecp = interp3_matrix(x1d, y1d, z1d, CP(:,1), CP(:,2), CP(:,3), porder, Band);

Eplot = interp3_matrix(x1d, y1d, z1d, Sphere.Location(:,1), Sphere.Location(:,2), Sphere.Location(:,3), porder, Band);

GCart = make3DImplicitGaussian(x1d, y1d, z1d, sigma, spacing, Band, numsigmas, LimitFarPoints);

G = GCart*Ecp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameter Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ScaleParameter = findScaleParamter(sigma, alpha, NumSteps, 'natural', '2d');

ScaleParameter = sqrt(0:NumSteps)'*sqrt(2)*sigma;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Signal = zeros(length(CPSignal), NumSteps);
% Signal(:,1) = CPSignal;
Signal = CPSignal;

Truth = zeros(Sphere.LocationCount, NumSteps);

SignalAtVertex = zeros(Sphere.LocationCount, NumSteps);
SignalAtVertex(:,1) = SignalOriginal;
% SignalAtVertex = SignalOriginal;

WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumSteps-1));
AbsErr = zeros(NumSteps,1);

if ShowPlot
    figure(1)
end

for i = 1 : NumSteps - 1
  
    SignalNew = G * Signal;

    SignalAtVertex(:,i+1) = Eplot * SignalNew;

    Truth= ExactSignal(ScaleParameter(i+1), Phi);
    AbsErr(i,1) = norm(Truth - SignalAtVertex(:,i+1), inf);

    if ShowPlot
        clf
        plot(Phi, SignalOriginal,'k.')
        hold on
        plot(Phi, SignalAtVertex(:,i),'r.')
        plot(Phi, Truth, 'gd')
    end
  
    Signal = SignalNew;
    
    waitbar(i/NumSteps, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumSteps-1));
end

waitbar(i/NumSteps, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)
if ShowPlot
    close(figure(1))
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% Andrew Rhodes
% ASEL
% February 2018

% Diffusion accuracy comparison for a signal on a sphere.


% close all
clear
clc


%% Additional Paths

addpath('~/GitProjects/matlab-utilities/')
addpath('~/Desktop/MLIDAR-1.0/MATLAB_Modules/')
addpath('~/AFOSR/Ashish/CS3 Code/')
addpath(genpath('~/GitProjects/pose/MATLAB_PointCloudDescriptors/OURCVFH/'))
addpath(genpath('~/Documents/Software/cp_matrices/'))
addpath(genpath('../'))
addpath('../src/')
addpath('../data/')
addpath(genpath('../models/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SpaceVec = [0.5, 0.25, 0.15, 0.1, 0.05];
% for MCRun = 1 : 5


NumberDivisions = 3;
alpha = 1;


FileLocation = '../models/Sphere/';
FileName = strcat('Icosphere',num2str(NumberDivisions),'.off');

porder = 4; 
dim = 3;
Lorder = 2;
spacing = 0.5;
bandwidth = 1.00001*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));

tauImplicit = spacing / 8;
MaxTauImplicit = 1/spacing;
NumStepsImplicit = 5000;%ceil(MaxTauImplicit / tauImplicit);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Sphere and Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sphere.Location,Sphere.Face] = icosphere(NumberDivisions);

Sphere.FaceCount = size(Sphere.Face, 1);
Sphere.LocationCount = size(Sphere.Location,1);
Sphere.FaceArea = findFaceAreas(Sphere.Location,Sphere.Face);



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

MinPoint = round(min(Sphere.Location) - bandwidth - spacing, 1);
MaxPoint = round(max(Sphere.Location) + bandwidth + spacing, 1);


[IJK,DIST,CP,XYZ,CPFACE] = tri2cp(Sphere.Face, Sphere.Location, spacing, MinPoint, porder, 1);

x1d = (MinPoint(1):spacing:MaxPoint(1))';
y1d = (MinPoint(2):spacing:MaxPoint(2))';
z1d = (MinPoint(3):spacing:MaxPoint(3))';

BandSearchSize = [length(y1d), length(x1d), length(z1d)];

Band = sub2ind(BandSearchSize, IJK(:,2), IJK(:,1), IJK(:,3));

FaceInterpolateWeights = interpBarycenterTriangle(Sphere, CP, CPFACE);

CPSignal = FaceInterpolateWeights * Sphere.Signal;

L = laplacian_3d_matrix(x1d,y1d,z1d, Lorder, Band);

Ecp = interp3_matrix(x1d, y1d, z1d, CP(:,1), CP(:,2), CP(:,3), porder, Band);

Eplot = interp3_matrix(x1d, y1d, z1d, Sphere.Location(:,1), Sphere.Location(:,2), Sphere.Location(:,3), porder, Band);


M = lapsharp(L, Ecp);

ItM = speye(size(M)) - alpha*tauImplicit * M;

I23tM = speye(size(M)) - (2/3)*alpha*tauImplicit * M;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameter Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ScaleParameter = zeros(NumStepsImplicit,1);

for i = 1 : NumStepsImplicit - 1
   
%     ScaleParameter(i+1,1) = sqrt(2*alpha*i*tauImplicit);
    ScaleParameter(i+1,1) = alpha*i*tauImplicit;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Signal = zeros(length(CPSignal), NumStepsImplicit);
Signal(:,1) = CPSignal;

SignalAtVertex = zeros(Sphere.LocationCount, NumStepsImplicit);
SignalAtVertex(:,1) = SignalOriginal;

WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsImplicit-1));

for i = 1 : NumStepsImplicit - 1
    
    if i == 1
        [Signal(:,i+1), flag] = gmres(ItM, Signal(:,i), [], 1e-10, 50);
    else
        [Signal(:,i+1), flag] = gmres(I23tM, (4/3)*Signal(:,i) - (1/3)*Signal(:,i-1), [], 1e-10, 50);
    end
    
    SignalAtVertex(:,i+1) = Eplot * Signal(:,i+1);
    
    if flag
        disp(flag)
    end
   
    waitbar(i/NumStepsImplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsImplicit-1));
    
end

waitbar(i/NumStepsImplicit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% AbsErr(MCRun,1) = norm( bsxfun(@minus, mean(SignalOriginal)/2 , SignalAtVertex(:,i)), inf);





% end



% figure
% loglog(1./SpaceVec, AbsErr)






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



















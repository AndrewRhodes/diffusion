% Andrew Rhodes
% ASEL
% March 2018

% Diffusion accuracy comparison for a signal on a sphere.


% close all
clear
clc

addprojectpaths % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCspacing = [0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025];


MCError = zeros(length(MCspacing), 2);
MCErrorAll = cell(length(MCspacing));




for MC = 1 : length(MCspacing)
    
    
    clearvars -except MCspacing MC MCError MCErrorAll
    spacing = MCspacing(MC)

    
    
alpha = 1;
porder = 3; 
dim = 2;
Lorder = 2;
% spacing = 0.01;
% sigma <= spacing
sigma = spacing;
numsigmas = 7;
LimitFarPoints = 0;

if spacing > sigma
    bandwidth = 1.00001*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));
else
    bandwidth = 1.00001*numsigmas*sigma*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));
end

NumberDivisions = 3;


FileLocation = '../models/Sphere/';
FileName = strcat('Icosphere',num2str(NumberDivisions),'.off');


% tauImplicit = spacing / 8;
% MaxTauImplicit = 1/spacing;
NumSteps = round(1/spacing); % 5000 ; %



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

FaceInterpolateWeights = interpBarycenterTriangle(Sphere, CP, CPFACE);

CPSignal = FaceInterpolateWeights * Sphere.Signal;

% L = laplacian_3d_matrix(x1d,y1d,z1d, Lorder, Band);

Ecp = interp3_matrix(x1d, y1d, z1d, CP(:,1), CP(:,2), CP(:,3), porder, Band);

Eplot = interp3_matrix(x1d, y1d, z1d, Sphere.Location(:,1), Sphere.Location(:,2), Sphere.Location(:,3), porder, Band);

G = make3DImplicitGaussian(x1d, y1d, z1d, sigma, spacing, Band, numsigmas, LimitFarPoints);
% G = make2DImplicitGaussian(x1d, y1d, sigma, spacing, band, numsigmas, LimitFarPoints);

% GE = diag(G) + (G - diag(G))*Ecp;
GE = G*Ecp;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameter Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ScaleParameter = findScaleParamter(sigma, alpha, NumSteps, 1, 2);

% ScaleParameter = zeros(NumSteps,1);

% ScaleParameter = zeros(NumSteps,1);
% 
% for i = 1 : NumSteps 
%    
%     ScaleParameter(i+1,1) = sqrt(i)*sigma;
% 
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Signal = zeros(length(CPSignal), NumSteps);
% Signal(:,1) = CPSignal;
Signal = CPSignal;

% SignalAtVertex = zeros(Circle.LocationCount, NumSteps);
% SignalAtVertex(:,1) = SignalOriginal;
SignalAtVertex = SignalOriginal;

WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumSteps-1));
AbsErr = zeros(NumSteps,1);

% figure(1)
for i = 1 : NumSteps - 1
  
%     Signal(:,i+1) = GE * Signal(:,i);
    SignalNew = GE * Signal;


%     SignalAtVertex(:,i+1) = Eplot * Signal(:,i+1);
    SignalAtVertex = Eplot * SignalNew;

    
%     clf
%     plot(Theta, SignalOriginal,'b')
%     hold on
%     plot(Theta, SignalAtVertex,'k')
    Truth = exp(-ScaleParameter(i,1)^2/2) .* SignalOriginal;
% % % %     Truth = exp(-ScaleParameter(i,1)^2/2) .* cos(Theta) + exp(-9*ScaleParameter(i,1)^2/2) .* sin(3*Theta);
%     plot(Theta, Truth,'r--')
    AbsErr(i+1,1) = norm(Truth - SignalAtVertex, inf);


    Signal = SignalNew;
   
    waitbar(i/NumSteps, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumSteps-1));
end

waitbar(i/NumSteps, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)
% close(figure(1))


MCError(MC, 1:2) = [NumSteps, AbsErr(NumSteps - 1)];
MCErrorAll{MC, 1} = [(1:NumSteps)', AbsErr];

end

close all




figure(2)
loglog(MCError(:,1), MCError(:,2),'bd-')
hold on
logx = [10,100];
logy = (10e-1).*logx.^(-2);
% logx = [100,1000];
% logy = (10e-1).*logx.^(-2);
loglog(logx, logy,'k-')




figure(3)
for i = 1 : length(MCErrorAll)
    loglog(MCErrorAll{i,1}(:,1), MCErrorAll{i,1}(:,2))
    hold on
end
% loglog(1:NumStepsImplicit, AbsErr)
xlabel('Iteration Number')
ylabel('Absolute Error')

logx = [10,1000];
logy = (10e-11).*logx.^(1);
loglog(logx, logy,'k-')







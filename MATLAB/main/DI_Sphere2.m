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

MaxDegreeL = 100;

porder = 5; % order of interpolation
dim = 3; % dimension
Lorder = 2; % Cartesian Laplace order
spacing = 0.1; % spacing of embedding grid
bandwidth = 1.002*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));

tauImplicit = spacing / 8; % time step
MaxTauImplicit = 1/spacing;
NumStepsImplicit = round(MaxTauImplicit); %ceil(MaxTauImplicit / tauImplicit);

ShowPlot = 1;

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

[xp, yp, zp] = sphere(65);
Sphere.Location(:,1) = xp(:);
Sphere.Location(:,2) = yp(:);
Sphere.Location(:,3) = zp(:);
Sphere.LocationCount = length(Sphere.Location);



[Sphere.Theta, Sphere.Phi, Sphere.Radius] = cart2sph(Sphere.Location(:,1), Sphere.Location(:,2), Sphere.Location(:,3));

SphericalHarmonic = makeRealSphericalHarmonic( MaxDegreeL, Sphere.Theta, Sphere.Phi );


% Define the signal
ExactSignal = @(sigma, SignalOriginal, MaxDegreeL) sum(cell2mat(cellfun(@times, num2cell( (exp(-(sigma^2/2).*(1:MaxDegreeL).*((1:MaxDegreeL)+1))')), SignalOriginal, 'UniformOutput', 0)),1);


SignalOriginal = ExactSignal(0, SphericalHarmonic, MaxDegreeL);
Sphere.Signal = SignalOriginal;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the Laplace Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MinPoint = min(Sphere.Location) - bandwidth - spacing;
MaxPoint = max(Sphere.Location) + bandwidth + spacing;

x1d = (MinPoint(1):spacing:MaxPoint(1))';
y1d = (MinPoint(2):spacing:MaxPoint(2))';
z1d = (MinPoint(3):spacing:MaxPoint(3))';

[GridX, GridY, GridZ] = meshgrid(x1d, y1d, z1d);
[CP(:,1), CP(:,2), CP(:,3), dist] = cpSphere(GridX(:), GridY(:), GridZ(:));


BandInit = find(abs(dist) <= bandwidth);

CPInit = CP(BandInit, :);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matric Construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ecp = interp3_matrix(x1d, y1d, z1d, CP(:,1), CP(:,2), CP(:,3), porder);
[EcpRow, EcpCol, EcpVal] = find(Ecp);
BandInner = unique(EcpCol);

L = laplacian_3d_matrix(x1d, y1d,z1d, Lorder, BandInner, BandInit);
[LRow, LCol, LVal] = find(L);
BandOuterTemp = unique(LCol);
BandOuter = BandInit( BandOuterTemp );

CPOut = CPInit(BandOuterTemp,:);

% Reform the L, Ecp matrices
Ecp = Ecp(BandOuterTemp, BandInner);
L = L(:, BandOuterTemp);



Eplot = interp3_matrix(x1d, y1d, z1d, Sphere.Location(:,1), Sphere.Location(:,2), Sphere.Location(:,3), porder, BandInner);


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


[CPTheta, CPPhi, CPRadius] = cart2sph(CPIn(:,1), CPIn(:,2), CPIn(:,3));

SphericalHarmonicCP = makeRealSphericalHarmonic( MaxDegreeL, CPTheta, CPPhi );

CPSignal = ExactSignal(0, SphericalHarmonicCP, MaxDegreeL);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M = lapsharp_unordered(L, Ecp, R);


ItM = speye(size(M)) - alpha*tauImplicit * M;

I23tM = speye(size(M)) - (2/3)*alpha*tauImplicit * M;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameter Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ScaleParameter = findScaleParamter(tauImplicit, alpha, NumStepsImplicit, 1, 3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Signal = zeros(length(CPSignal), NumStepsImplicit);
Signal(:,1) = CPSignal;

SignalAtVertex = zeros(Sphere.LocationCount, NumStepsImplicit);
SignalAtVertex(:,1) = SignalOriginal;

AbsErr = zeros(NumStepsImplicit,1);


WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsImplicit-1));
if ShowPlot
    figure(1)
end

for i = 1 : NumStepsImplicit - 1
    
%     if i == 1
%     Signal(:,i+1) = ItM \ Signal(:,i);
        [Signal(:,i+1), flag] = gmres(ItM, Signal(:,i), [], 1e-10, 100);
%     else
%         [Signal(:,i+1), flag] = gmres(I23tM, (4/3)*Signal(:,i) - (1/3)*Signal(:,i-1), [], 1e-10, 100);
%     end
    
    % Interpolate back to explicit surface
    SignalAtVertex(:,i+1) = Eplot * Signal(:,i+1);
    
    if flag
        disp(flag)
    end
    
    % Calculate Truth and Error
    Truth = ExactSignal(ScaleParameter(i), SphericalHarmonic, MaxDegreeL);
    AbsErr(i,1) = norm(Truth - SignalAtVertex(:,i), inf);
    
    
    if ShowPlot       
        clf
        plot(Sphere.Theta, SignalOriginal,'ko')
        hold on
        plot(Sphere.Theta, Truth, 'gd')
        plot(Sphere.Theta, SignalAtVertex(:,i),'r.')
        legend('Original', 'Truth at i', 'Diffused at i')

    end
   
    waitbar(i/NumStepsImplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsImplicit-1));
    
end

waitbar(i/NumStepsImplicit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)
if ShowPlot
    close(figure(1))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure
loglog(1:NumStepsImplicit, AbsErr)
xlabel('Iteration Number')
ylabel('Relative Error')

figure
plot(1:NumStepsImplicit, AbsErr)
xlabel('Iteration Number')
ylabel('Relative Error')



















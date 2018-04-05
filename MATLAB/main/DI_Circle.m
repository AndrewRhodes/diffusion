% Andrew Rhodes
% ASEL
% March 2018

% Diffusion accuracy comparison for a signal on a circle.


close all
clear
clc

addpath('../src/')
ProjectRoot = setupprojectpaths % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


alpha = 1;
porder = 4;
dim = 2;
Lorder = 2;
spacing = 0.01;
bandwidth = 1.002*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));

tauImplicit = spacing / 8;
MaxTauImplicit = 1/spacing;
NumStepsImplicit = round(MaxTauImplicit); 

% ExactSignal = @(sigma, theta) exp(-sigma^2/2)*cos(theta) + exp(-9*sigma^2/2)*sin(3*theta);
ExactSignal = @(sigma, theta) exp(-sigma^2/2)*cos(theta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Circle and Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Theta = linspace(0, 2*pi, 2000)';

Radius = ones(size(Theta));
[xp,yp] = pol2cart(Theta, Radius);
Circle.Location(:,1) = xp(:);
Circle.Location(:,2) = yp(:);
Circle.LocationCount = length(Circle.Location);


SignalOriginal = ExactSignal(0, Theta); 
Circle.Signal = SignalOriginal;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Grid and Inner/Outer Bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MinPoint = min(Circle.Location) - bandwidth - spacing;
MaxPoint = max(Circle.Location) + bandwidth + spacing;

x1d = (MinPoint(1):spacing:MaxPoint(1))';
y1d = (MinPoint(2):spacing:MaxPoint(2))';
        
        
% x1d = (-2:spacing:2)';
% y1d = x1d;%(-1.4:spacing:1.4)';

[GridX, GridY] = meshgrid(x1d, y1d);

[CP(:,1), CP(:,2), dist] = cpCircle(GridX(:), GridY(:));


BandInit = find(abs(dist) <= bandwidth);

CPInit = CP(BandInit,:);

GridXBandInit = GridX(BandInit);
GridYBandInit = GridY(BandInit);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matric Construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [L, Ecp, R, BandInner, BandOuter, BandInner2, BandOuter2] = ...
%     ops_and_bands2d(x1d, y1d, GridXBandInit, GridYBandInit, ...
%                     CPInit(:,1), CPInit(:,2), BandInit, porder, Lorder);
% 
% 
% CPIn = CPInit(BandInner,:);      
                
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


Eplot = interp2_matrix(x1d, y1d, Circle.Location(:,1), Circle.Location(:,2), porder, BandInner);
% Eplot = interp2_matrix(x1d, y1d, Circle.Location(:,1), Circle.Location(:,2), porder, BandInner2);

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
GridXIn = R*GridXOut;
GridYIn = R*GridYOut;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the Signal and Plot Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[CPTheta, CPr] = cart2pol(CPIn(:,1),CPIn(:,2));

CPSignal = ExactSignal(0, CPTheta); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the Laplace-Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M = lapsharp_unordered(L, Ecp, R);

ItM = speye(size(M)) - alpha*tauImplicit * M;

I23tM = speye(size(M)) - (2/3)*alpha*tauImplicit * M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameter Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ScaleParameter = findScaleParamter(tauImplicit, alpha, NumStepsImplicit, 1, 3);


ScaleParameter = zeros(NumStepsImplicit,1);

for i = 1 : NumStepsImplicit - 1

%     ScaleParameter(i+1,1) = sqrt(2*i*alpha*tauImplicit);
%     ScaleParameter(i+1,1) = alpha*i*tauImplicit;
    ScaleParameter(i+1,1) = sqrt(2*i*alpha*tauImplicit);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Signal = zeros(length(CPSignal), NumStepsImplicit);
Signal(:,1) = CPSignal;

SignalAtVertex = zeros(Circle.LocationCount, NumStepsImplicit);
SignalAtVertex(:,1) = SignalOriginal;

WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsImplicit-1));
AbsErr = zeros(NumStepsImplicit,1);
figure(1)
for i = 1 : NumStepsImplicit - 1

    if i == 1
        Signal(:,i+1) = ItM \ Signal(:,i);
%         [Signal(:,i+1), flag] = gmres(ItM, Signal(:,i), [], 1e-10, 100);
    else
        Signal(:,i+1) = I23tM \ ((4/3)*Signal(:,i) - (1/3)*Signal(:,i-1));
% %         [Signal(:,i+1), flag] = gmres(I23tM, (4/3)*Signal(:,i) - (1/3)*Signal(:,i-1), 2, 1e-10, 50);
    end

    SignalAtVertex(:,i+1) = Eplot * Signal(:,i+1);

    clf
    plot(Theta, SignalOriginal,'b')
    hold on
    plot(Theta, SignalAtVertex(:,i+1),'k')

    Truth = ExactSignal(ScaleParameter(i), Theta); 
    plot(Theta, Truth,'r--')
    AbsErr(i+1,1) = norm(Truth - SignalAtVertex(:,i+1), inf);

    if flag
        disp(flag)
    end

    waitbar(i/NumStepsImplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsImplicit-1));
    pause(0.2)
end

waitbar(i/NumStepsImplicit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)
close(figure(1))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
loglog(1:NumStepsImplicit,AbsErr,'k-')
xlabel('Iteration Number')
ylabel('Absolute Error')





% SignalAtVertex(abs(SignalAtVertex)<eps) = 0;
%
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








%
% figure
% loglog(MCError(:,1), MCError(:,2),'bd-')
% hold on
% logx = [5,100];
% logy = (10e-4).*logx.^(-2);
% loglog(logx, logy,'k-')
%
%
%
%
% figure
% for i = 1 : length(MCErrorAll)
%     loglog(MCErrorAll{i,1}(:,1), MCErrorAll{i,1}(:,2))
%     hold on
% end
% % loglog(1:NumStepsImplicit, AbsErr)
% xlabel('Iteration Number')
% ylabel('Absolute Error')

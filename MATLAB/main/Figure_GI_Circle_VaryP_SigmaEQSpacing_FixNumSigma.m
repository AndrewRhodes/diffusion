% Andrew Rhodes
% ASEL
% March 2018

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

MCspacing = [0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001];

MCporder = [1, 2, 3, 4, 5, 6, 7];

MCError = zeros(length(MCspacing), 2, length(MCporder));
MCErrorAll = cell(length(MCspacing), length(MCporder));


for MCp = 1 : length(MCporder)

for MCs = 1 : length(MCspacing)
    
    
    clearvars -except MCspacing MCs MCp MCError MCErrorAll MCporder
    spacing = MCspacing(MCs)
    porder = MCporder(MCp)
    
    
alpha = 1;
% porder = 9; 
dim = 2;
Lorder = 2;
% spacing = 0.01;
% sigma <= spacing
sigma = spacing;
numsigmas = 5;
LimitFarPoints = 1;

if spacing > sigma
    bandwidth = 1.00001*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));
else
    bandwidth = 1.00001*numsigmas*sigma*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));
end


% tauImplicit = spacing / 8;
% MaxTauImplicit = 1/spacing;
NumSteps = round(1/spacing); % 5000 ; %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Circle and Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Theta = linspace(0,2*pi,200)';

Radius = ones(size(Theta));
[xp,yp] = pol2cart(Theta, Radius);
Circle.Location(:,1) = xp(:);
Circle.Location(:,2) = yp(:);
Circle.LocationCount = length(Circle.Location);


SignalOriginal = cos(Theta) + sin(3*Theta);

Circle.Signal = SignalOriginal;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the Laplace Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MinPoint = round(min(Circle.Location(:,1)) - bandwidth - ceil(numsigmas * sigma), 1);
MaxPoint = round(max(Circle.Location(:,1)) + bandwidth + ceil(numsigmas * sigma), 1);

% MinPoint = round(min(Circle.Location(:,1)) - bandwidth, 1);
% MaxPoint = round(max(Circle.Location(:,1)) + bandwidth, 1);


x1d = (MinPoint:spacing:MaxPoint)';
y1d = x1d;

[GridX, GridY] = meshgrid(x1d, y1d);

[CP(:,1), CP(:,2), dist] = cpCircle(GridX(:), GridY(:));


% outer_band = find(abs(dist) <= 2*bandwidth);
band = find(abs(dist) <= bandwidth);

CP = CP(band,:);

GridXBand = GridX(band); 
GridYBand = GridY(band);


[CPTheta, CPr] = cart2pol(GridXBand,GridYBand);


CPSignal = cos(CPTheta) + sin(3*CPTheta);


Ecp = interp2_matrix(y1d, x1d, CP(:,2), CP(:,1), porder, band);
% Ecp = interp2_matrix(x1d, y1d, CP(:,1), CP(:,2), porder, band);
% Ecp = Ecp(outer_band,inner_band);
% size(Ecp)



Eplot = interp2_matrix(x1d, y1d, Circle.Location(:,1), Circle.Location(:,2), porder, band);

% G = make3DImplicitGaussian(y1d, x1d, z1d, sigma, spacing, band, 4, 1);
G = make2DImplicitGaussian(x1d, y1d, sigma, spacing, band, numsigmas, LimitFarPoints);


% GE = diag(G) + (G - diag(G))*Ecp;
GE = G*Ecp;


% Gaussian Method
% MinPoint = round(min(PointCloud.Location) - bandwidth - 2*spacing, 1);
% MaxPoint = round(max(PointCloud.Location) + bandwidth + 2*spacing, 1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameter Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ScaleParameter = findScaleParamter(sigma, alpha, NumSteps, 1, 2);

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
%     Truth = exp(-ScaleParameter(i,1)^2/2) .* SignalOriginal;
    Truth = exp(-ScaleParameter(i,1)^2/2) .* cos(Theta) + exp(-9*ScaleParameter(i,1)^2/2) .* sin(3*Theta);
%     plot(Theta, Truth,'r--')
    AbsErr(i+1,1) = norm(Truth - SignalAtVertex, inf);


    Signal = SignalNew;
   
    waitbar(i/NumSteps, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumSteps-1));
end

waitbar(i/NumSteps, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)
% close(figure(1))


MCError(MCs, 1:2, MCp) = [NumSteps, AbsErr(NumSteps - 1)];
MCErrorAll{MCs, MCp} = [(1:NumSteps)', AbsErr];

end
end

close all




% figure(2)
% loglog(MCError(:,1), MCError(:,2),'bd-')
% hold on
% logx = [100,10000];
% logy = (10e-1).*logx.^(-2);
% loglog(logx, logy,'k-')
% 
% 
% 
% 
% figure(3)
% for i = 1 : length(MCErrorAll)
%     loglog(MCErrorAll{i,1}(:,1), MCErrorAll{i,1}(:,2))
%     hold on
% end
% % loglog(1:NumStepsImplicit, AbsErr)
% xlabel('Iteration Number')
% ylabel('Absolute Error')
% 
% logx = [10,1000];
% logy = (10e-11).*logx.^(1);
% loglog(logx, logy,'k-')




% MCErrorAll{MCs, MCp} = [(1:NumSteps)', AbsErr];

colors = ['k','b','r','c','m','g','k'];
shapes = ['.','^','.','*','+','o','s'];
lin  = {':','-','-','-','-.',':','--'};
msizes = [8, 8, 8, 8, 10, 12, 15];

Points = zeros(size(MCErrorAll,1),2, size(MCErrorAll,2));
for j = 1 : size(MCErrorAll,2)
    for i = 1 : size(MCErrorAll,1)
        Points(i,1:2,j) = [MCErrorAll{i,j}(end,1), MCErrorAll{i,j}(end,2)];
    end
end


figure('units','normalized','outerposition',[0 0 1 1])
for i = 1 : size(Points,3)
    loglog( Points(:,1,i), Points(:,2,i), strcat(colors(i), shapes(i), lin{i}), 'linewidth', 3, 'markersize', msizes(i) )
    hold on
end

logx = [10,100];
logy = (10e-2).*logx.^(-2);
loglog(logx, logy,'k-','linewidth',3)


logx = [100,1000];
logy = (10*10^(-0.5)).*logx.^(-1);
loglog(logx, logy,'k-','linewidth',3)


xl=get(gca,'XLim');
yl=get(gca,'YLim');
text1 = text(10,10*10^(-4.5),'$2^{nd}$ Order','Interpreter','latex');
set(text1,'Rotation',-32)
set(text1,'FontSize',50)


text2 = text(200,10*10^(-2.5),'$1^{st}$ Order','Interpreter','latex');
set(text2,'Rotation',-16)
set(text2,'FontSize',50)


xlabel('N')
ylabel('$\| $error$ \|_{\infty}$','Interpreter','latex')
ax = gca;
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;


xlim([10*10^(-0.75) 10*10^(2.5)])

xticks([10e0 10e1 10e2 10e3])
yticks([10e-6, 10e-5, 10e-4, 10e-3, 10e-2, 10e-1, 10e0])
yticklabels({'10e-6','10e-5', '10e-4', '10e-3', '10e-2', '10e-1', '10e0'})

hleg = legend('p=1','p=2','p=3','p=4','p=5','p=6','p=7');
set(hleg, 'position', [0.85 0.72 0.05 0.2])













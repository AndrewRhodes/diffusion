% Andrew Rhodes
% ASEL
% March 2018

% Diffusion accuracy comparison for a signal on a sphere.


close all
clear
clc

ProjectRoot = addprojectpaths % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCspacing = [1/3, 1/9, 1/22, 1/30, 1/90, 1/220, 1/300];
% MCspacing = [0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025];%, 0.001];
MCporder = [1, 2, 3, 4, 5, 6, 7];

MCError = zeros(length(MCspacing), 2, length(MCporder));
MCErrorAll = cell(length(MCspacing), length(MCporder));


for MCp = 1 : length(MCporder)

for MCs = 1 : length(MCspacing)
    
    
    clearvars -except MCspacing MCporder MCs MCp MCError MCErrorAll 
    spacing = MCspacing(MCs)
    porder = MCporder(MCp)
    
    
alpha = 1;

% porder = 5; 
dim = 2;
Lorder = 2;
% spacing = 0.01;
bandwidth = 1.00001*spacing*sqrt((dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2));

tauImplicit = spacing / 8;
MaxTauImplicit = 1/spacing;
NumStepsImplicit = round(MaxTauImplicit); % 3000;%100*MaxTauImplicit; %ceil(MaxTauImplicit / tauImplicit);

% NumStepsImplicit = 500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Circle and Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Theta = linspace(0,2*pi,1000)';

Radius = ones(size(Theta));
[xp,yp] = pol2cart(Theta, Radius);
Circle.Location(:,1) = xp(:);
Circle.Location(:,2) = yp(:);
% Circle.Location(:,3) = zeros(size(yp(:)));
Circle.LocationCount = length(Circle.Location);


SignalOriginal = cos(Theta) + sin(3*Theta);
% SignalOriginal = cos(Theta);


Circle.Signal = SignalOriginal;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the Laplace Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MinPoint = (min(Circle.Location) - bandwidth - spacing);
MaxPoint = (max(Circle.Location) + bandwidth + spacing);

x1d = (MinPoint(1):spacing:MaxPoint(1))';
y1d = (MinPoint(2):spacing:MaxPoint(2))';
% z1d = (MinPoint(3):spacing:MaxPoint(3))';

[GridX, GridY] = meshgrid(x1d, y1d);
% [GridX, GridY, GridZ] = meshgrid(x1d, y1d, z1d);

[CP(:,1), CP(:,2), dist] = cpCircle(GridX(:), GridY(:));


% outer_band = find(abs(dist) <= 2*bandwidth);
band = find(abs(dist) <= bandwidth);

CP = CP(band,:);

GridXBand = GridX(band); 
GridYBand = GridY(band);


[CPTheta, CPr] = cart2pol(CP(:,1),CP(:,2));


CPSignal = cos(CPTheta)+ sin(3*CPTheta);
% CPSignal = cos(CPTheta);

Ecp = interp2_matrix(x1d, y1d, CP(:,1), CP(:,2), porder, band);
% Ecp = Ecp(outer_band,inner_band);
% size(Ecp)


L = laplacian_2d_matrix(x1d, y1d, Lorder, band);
% size(L)


Eplot = interp2_matrix(x1d, y1d, Circle.Location(:,1), Circle.Location(:,2), porder, band);


% spy(L(inner_band,:))
% innerouter= zeros(length(inner_band),1);
% 
% for i = 1 : length(inner_band)
%    [a,b] = find(outer_band == inner_band(i));
%    if ~isempty(a)
%        innerouter(i) = 1;
%    end
%        
%    
% end


M = lapsharp(L, Ecp);

ItM = speye(size(M)) - alpha*tauImplicit * M;

I23tM = speye(size(M)) - (2/3)*alpha*tauImplicit * M;

I611M = speye(size(M)) - (6/11)*alpha*tauImplicit * M;

I1225M = speye(size(M)) - (12/25)*alpha*tauImplicit * M;

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

NumIter = min(length(band),100);

% figure(1)
for i = 1 : NumStepsImplicit - 1
  
%     if i == 1
        Signal(:,i+1) = ItM \ Signal(:,i);
%         [Signal(:,i+1), flag] = gmres(ItM, Signal(:,i), [], 1e-10, NumIter);
%     elseif i == 2
% %         Signal(:,i+1) = I23tM \ ((4/3)*Signal(:,i) - (1/3)*Signal(:,i-1));
%         [Signal(:,i+1), flag] = gmres(I23tM, (4/3)*Signal(:,i) - (1/3)*Signal(:,i-1), [], 1e-10, NumIter);
%     elseif i == 3
% %         Signal(:,i+1) = I611M \ ((4/3)*Signal(:,i) - (1/3)*Signal(:,i-1));
%         [Signal(:,i+1), flag] = gmres(I611M, (18/11)*Signal(:,i) - (9/11)*Signal(:,i-1) + (2/11)*Signal(:,i-2), [], 1e-10, NumIter);    
%     else
% %         Signal(:,i+1) = I1225M \ ((4/3)*Signal(:,i) - (1/3)*Signal(:,i-1));
%         [Signal(:,i+1), flag] = gmres(I1225M, (48/25)*Signal(:,i)-(36/25)*Signal(:,i-1) + (16/25)*Signal(:,i-2) - (3/25)*Signal(:,i-3), [], 1e-10, NumIter);    
%     end
    
    SignalAtVertex(:,i+1) = Eplot * Signal(:,i+1);
    
    clf
    plot(Theta, SignalOriginal,'b')
    hold on
    plot(Theta, SignalAtVertex(:,i+1),'k')
% % % %     Truth = exp(-ScaleParameter(i+1)^2/2) .* SignalOriginal;
    Truth = exp(-ScaleParameter(i+1)^2/2) .* cos(Theta) + exp(-9*ScaleParameter(i+1)^2/2) .* sin(3*Theta);
    plot(Theta, Truth,'r--')
    AbsErr(i+1,1) = norm(Truth - SignalAtVertex(:,i+1), inf);
    
    if flag
        disp(flag)
    end
   
    waitbar(i/NumStepsImplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsImplicit-1));
%     pause(0.2)
end

waitbar(i/NumStepsImplicit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)
% close(figure(1))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SignalAtVertex(abs(SignalAtVertex)<eps) = 0;

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

% MCError(MC,1:2) = [NumStepsImplicit, AbsErr(NumStepsImplicit - 1)];
% MCErrorAll{MC,1} = [(1:NumStepsImplicit)', AbsErr];

MCError(MCs, 1:2, MCp) = [NumStepsImplicit, AbsErr(NumStepsImplicit - 1)];
MCErrorAll{MCs, MCp} = [(1:NumStepsImplicit)', AbsErr];

end
end


close all


colors = ['k','k','b','r','r','c','m','y','g','k'];
shapes = {'','','^','h','p','>','*','+','o','s'};
lin  = {':','-.','-','-','-','-.','--','-.',':','--'};
msizes = [8, 8, 8, 8, 8, 8, 10, 12, 14, 17];


Points = zeros(size(MCErrorAll,1), 2, size(MCErrorAll,2));
for j = 1 : size(MCErrorAll,2)
    for i = 1 : size(MCErrorAll,1)
        Points(i,1:2,j) = [MCErrorAll{i,j}(end,1), MCErrorAll{i,j}(end,2)];
    end
end


figure('units','normalized','outerposition',[0 0 1 1])
for i = 1 : size(Points,3)
    loglog( Points(:,1,i), Points(:,2,i), strcat(colors(i), shapes{i}, lin{i}), 'linewidth', 3, 'markersize', msizes(i) )
    pause
    hold on
end


logx = [10,100];
logy = (10e-2).*logx.^(-2);
loglog(logx, logy,'k-','linewidth',3)



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











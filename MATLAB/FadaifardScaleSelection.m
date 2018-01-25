



close all
clear
clc


NumSteps = 50;

% Find the scale parameters
maxsample = 10; 
ws = 0 : 0.05 : maxsample;
NumSample = length(ws);
ScaleParameterSpatial1 = zeros(NumSteps,1);
ScaleParameterSpatial2 = zeros(NumSteps,1);
ScaleParameterFrequency1 = zeros(NumSteps,1);
ScaleParameterFrequency2 = zeros(NumSteps,1);
CutoffFrequency1 = zeros(NumSteps,1);
CutoffFrequency2 = zeros(NumSteps,1);
cut = sqrt(log(2));
db3 = 1/sqrt(2);
ws2 = (ws.^2)';
tau = 1;
hEdge = 3;
H1 = zeros(NumSample,NumSteps);
H1(:,1) =  ones(NumSample,1) ;
H2 = zeros(NumSample,NumSteps);
H2(:,1) =  ones(NumSample,1) ;
h1 = 1 ./ ( ones(NumSample,1) + tau / hEdge^2 * ws2 );
h2 = 1 ./ ( ones(NumSample,1) + tau * ws2 );

for i = 1 : NumSteps - 1
    
    % Transfer function
    % 1st order
    H1(:,i+1) = H1(:,i) .* h1;
    % Find the frequency at the cutoff values
    CutoffFrequency1(i+1,1) = interp1(H1(:,i+1), ws, db3);
    % Change cutoff frequency to scale parameter
    ScaleParameterFrequency1(i,1) = CutoffFrequency1(i+1,1) / cut;   
    ScaleParameterSpatial1(i,1) = hEdge / ScaleParameterFrequency1(i,1);
    
    
    
    % Transfer function
    % 1st order
    H2(:,i+1) = H2(:,i) .* h2;
    % Find the frequency at the cutoff values
    CutoffFrequency2(i+1,1) = interp1(H2(:,i+1), ws, db3);
    % Change cutoff frequency to scale parameter
    ScaleParameterFrequency2(i,1) = CutoffFrequency2(i+1,1) / cut;   
    ScaleParameterSpatial2(i,1) = 1 / ScaleParameterFrequency2(i,1);
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fadiafard Method (with a constrained search range)





% fs_spatial = sqrt(tau) ./ (2*pi*fs_frequency);

j = 3;
figure
plot(ws, H2(:,j),'k.','markersize',12')
hold on
plot([CutoffFrequency2(j) CutoffFrequency2(j)],[0 db3],'r-','linewidth',5)
plot([0 CutoffFrequency2(j)],[db3 db3],'r-','linewidth',5)
plot(CutoffFrequency2(j), db3,'r.','markersize',10')
% plot(ws, exp(-0.5*ws2/ScaleParameterFrequency2(j,1)^2),'b-')
plot(ws, exp(-0.5*ws2 * (2*j*tau)),'b-','linewidth',5)

maxsamplef = [0.4082, 4, 8];

for n = maxsamplef
    
    wsf = (0:0.05:n)';
    wsf2 = wsf.^2;
    wsf4 = wsf2.^2;
    fs_frequency = zeros(NumSteps,1);
    
    for i = 1 : NumSteps
        fs_frequency(i,1) = sqrt( 0.5*sum(wsf4') ./ sum( (wsf2').*sum( log(1 + ones(i,1) * (wsf2')) ,1) ,2) );
    end
    
    plot(ws, exp(-0.5*ws2 / fs_frequency(j)^2), 'm-','linewidth',5)
    
    
end

xlim([0 3])
ylabel('Response')
xlabel('Spatial Frequency')
ax = gca;
ax.XAxis.FontSize = 55;
ax.YAxis.FontSize = 55;
tb1 = annotation('arrow',[0.55 0.35],[0.7 0.35]);
tb1.LineWidth = 5;
tb1.HeadWidth = 25;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot in log-space

% figure
% j = 3;
% plot(ws,log(H1(:,j)),'k.')
% hold on
% plot(ws, log(exp(-0.5*ws2 / ScaleParameterFrequency1(j,1)^2)),'b-')
% plot(ws, log(exp(-0.5*ws2 * (2*j*tau))), 'g--')
% xlim([0 maxsample])
% ylim([-70 0])
% legend('1st Order Implicit Euler TF','Least Squares Gaussian','Cutoff Frequency Gaussian')
% ylabel('Intensity')
% xlabel('Spatial Frequency')
% 
% 
% 
% j = 3;
% figure
% plot(ws,log(H2(:,j)),'k:')
% hold on
% plot(ws, log(exp(-0.5*ws2 / ScaleParameterFrequency2(j,1)^2)),'b-')
% plot(ws, log(exp(-0.5*ws2 * (2*j*tau))),'g--')
% axis([0 4 -20 0])






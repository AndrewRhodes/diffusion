% Andrew Rhodes
% WVU
% Jan. 2018
% Compare the scale estimation technique by Fadaifard to our technique



close all
clear
clc


NumSteps = 50;
tau = 1;

%% Method of cutoff frequency 

maxsample = 10; 
ws = 0 : 0.01 : maxsample;
NumSample = length(ws);
ScaleParameterSpatial = zeros(NumSteps,1);
ScaleParameterFrequency = zeros(NumSteps,1);
cut = sqrt(log(2));
db3 = 1/sqrt(2);
ws2 = (ws.^2)';

H = zeros(NumSample, NumSteps);
H(:,1) = ones(NumSample,1) ;
h = 1 ./ ( ones(NumSample,1) + tau * ws2 );
CutoffFrequency = zeros(NumSteps,1);


for i = 1 : NumSteps - 1
    
    % Transfer function
    % 1st order
    H(:,i+1) = H(:,i) .* h;
    % Find the frequency at the cutoff values
    [uH, indH] = unique(H(:,i+1));
    CutoffFrequency(i+1,1) = interp1(uH, ws(indH), db3);
    % Change cutoff frequency to scale parameter
    ScaleParameterFrequency(i+1,1) = CutoffFrequency(i+1,1) / cut;   
    ScaleParameterSpatial(i+1,1) =  1 / ScaleParameterFrequency(i+1,1);

end


%% Method of natural scale growth 


ScaleParameterSpatialNatural = zeros(NumSteps,1);

for i = 1 : NumSteps - 1
   
    ScaleParameterSpatialNatural(i+1,1) = sqrt(2*i*tau);

end





j = 2;
figure
plot(ws, H(:,j),'k--','linewidth',8)
hold on
plot([CutoffFrequency(j) CutoffFrequency(j)],[0 db3],'r-','linewidth',5)
plot([0 CutoffFrequency(j)],[db3 db3],'r-','linewidth',5)
plot(CutoffFrequency(j), db3,'r.','markersize',14)
plot(ws, exp(-0.5*ws2/ScaleParameterFrequency(j+1,1)^2),'b-','linewidth',5)
plot(ws, exp(-0.5*ws2 * (2*j*tau)),'g--','linewidth',5)
xlim([0 3])
xlabel('Spatial Frequency')
ylabel('Response')
ax = gca;
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;



%% Method of Faidarfard 

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






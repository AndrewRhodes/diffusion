% Andrew Rhodes
% WVU
% Dec. 2018
% Compare the scale estimation technique by Fadaifard to our technique
% For journal paper. Fig. 2

close all
clear
clc

NumSteps = 500;
tau = 1;
alpha = 1;

MaxFreq = 2;
StepFreq = 0.0001;


%% Natural Scale Growth Method


% Next two methods are exactly the same.
% t_n = cumsum([0;ones(NumSteps-1,1)]);
% sigma_n = sqrt(2*t_n);

% ScaleNatural = sqrt(2*alpha*tau*(0:NumSteps-1)');
ScaleNatural = sqrt(2*alpha*cumsum([0; ones(NumSteps-1,1)]).*(0:NumSteps-1)');

ScaleNaturalFreq = 1./ScaleNatural;


%% Cutoff Frequency Method


CutoffFrequency = zeros(NumSteps,1);
ScaleCutoffFreq = zeros(NumSteps,1);
ScaleCutoff = zeros(NumSteps,1);

db3 = 1/sqrt(2);
cut = sqrt(2*log(sqrt(2)));
ws = 0 : StepFreq : 5;
ws2 = (ws.^2)';
NumSample = length(ws);

H = zeros(NumSample, NumSteps);
H(:,1) = ones(NumSample, 1);
h = 1 ./ ( ones(NumSample,1) + alpha * tau * ws2);


for i = 1 : NumSteps - 1
    
    h = 1 ./ ( ones(NumSample,1) + alpha * i * ws2); 
    
    H(:,i+1) = H(:,i) .* h;
    % Find the frequency at the cutoff values
    [uH, indH] = unique(H(:,i+1));
    CutoffFrequency(i+1,1) = interp1(uH, ws(indH), db3);
    
    ScaleCutoffFreq(i+1,1) = abs(CutoffFrequency(i+1,1) / cut);   
    ScaleCutoff(i+1,1) =  1 / ScaleCutoffFreq(i+1,1);    
           
end




%% Quadratic Fitting Method (Fadaifard)


ScaleQuadFreq = zeros(NumSteps,1);
ScaleQuad1 = zeros(NumSteps,1);

% Test with omega_J = 2
ws = (-MaxFreq : StepFreq : MaxFreq)';
ws2 = ws.^2;
ws4 = ws2.^2;
sumws4 = sum(ws4);
delta = 1.1;

for i = 1 : NumSteps - 1 
        
    %% Fadaifard Method with exponential growing time step
    
%     ScaleQuadFreq(i+1,1) = sqrt( 0.5*sumws4 ...
%         ./ sum( (ws2').*sum( log(1 + tau * delta.^(1:i)' * (ws2')), 1) ,2) );

    
    %% Current Discussion of holding constant tau = 1
    
%     ScaleQuadFreq(i+1,1) = sqrt( 0.5*sumws4 ...
%         ./ sum( (ws2').*sum( log(1 + alpha * cumsum([0; ones(i-1,1)]) * (ws2')), 1) ,2) );
    
    ScaleQuadFreq(i+1,1) = sqrt( 0.5*sumws4 ...
        ./ sum( (ws2').*sum( log(1 + alpha * ones(i,1) * (ws2')), 1) ,2) );
    
    
    %% Convert to Spatial Domain
    
    ScaleQuad1(i+1,1) = 1 / (ScaleQuadFreq(i+1,1));
    
end



ScaleQuad2 = zeros(NumSteps,1);
% Test with omega_J = 4
ws = (-4 : StepFreq : 4)';
ws2 = ws.^2;
ws4 = ws2.^2;
sumws4 = sum(ws4);

for i = 1 : NumSteps - 1 
        
    %% Fadaifard Method with exponential growing time step
    
%     ScaleQuadFreq(i+1,1) = sqrt( 0.5*sumws4 ...
%         ./ sum( (ws2').*sum( log(1 + tau * delta.^(1:i)' * (ws2')), 1) ,2) );
    

    %% Current Discussion of holding constant tau = 1
%     ScaleQuadFreq(i+1,1) = sqrt( 0.5*sumws4 ...
%         ./ sum( (ws2').*sum( log(1 + alpha * cumsum([0; ones(i-1,1)]) * (ws2')), 1) ,2) );
%     
  ScaleQuadFreq(i+1,1) = sqrt( 0.5*sumws4 ...
    ./ sum( (ws2').*sum( log(1 + alpha * ones(i,1) * (ws2')), 1) ,2) );

    
    %% Convert to Spatial Domain
    
    ScaleQuad2(i+1,1) = 1 / (ScaleQuadFreq(i+1,1));
end




ScaleQuad3 = zeros(NumSteps,1);
% Test with omega_J = 5
ws = (-5 : StepFreq : 5)';
ws2 = ws.^2;
ws4 = ws2.^2;
sumws4 = sum(ws4);

for i = 1 : NumSteps - 1 
        

    %% Fadaifard Method with exponential growing time step
    
%     ScaleQuadFreq(i+1,1) = sqrt( 0.5*sumws4 ...
%         ./ sum( (ws2').*sum( log(1 + tau * delta.^(1:i)' * (ws2')), 1) ,2) );
    

    %% Current Discussion of holding constant tau = 1
      
%     ScaleQuadFreq(i+1,1) = sqrt( 0.5*sumws4 ...
%         ./ sum( (ws2').*sum( log(1 + alpha * cumsum([0; ones(i-1,1)]) * (ws2')), 1) ,2) );
    
    ScaleQuadFreq(i+1,1) = sqrt( 0.5*sumws4 ...
        ./ sum( (ws2').*sum( log(1 + alpha * ones(i,1) * (ws2')), 1) ,2) );
    
    
    %% Convert to Spatial Domain
    
    ScaleQuad3(i+1,1) = 1 / (ScaleQuadFreq(i+1,1));
end





%% Errors compared to natural scale growth


ErrorCutoff =  abs(ScaleNatural - ScaleCutoff) ./ ScaleNatural;
ErrorQuad1 =  abs(ScaleNatural - ScaleQuad1/(1-0.3030) ) ./ ScaleNatural;
ErrorQuad2 =  abs(ScaleNatural- ScaleQuad2/(1-0.5159) ) ./ ScaleNatural;
ErrorQuad3 =  abs(ScaleNatural- ScaleQuad3/(1-0.5797) ) ./ ScaleNatural;


figure
semilogy(ScaleNatural, ErrorCutoff, 'k', 'Linewidth', 4)
hold on
semilogy(ScaleNatural, ErrorQuad1, 'b--', 'Linewidth', 4)
semilogy(ScaleNatural, ErrorQuad2, 'r-.', 'Linewidth', 4)
semilogy(ScaleNatural, ErrorQuad3, 'g:', 'Linewidth', 4)
ylabel('Rel. Error', 'Fontsize', 50)
ax = gca;
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;
ax.YTick = [10^-4, 10^-2, 10^0];
xlabel('$\sigma_n$', 'Interpreter', 'Latex', 'Fontsize', 80)
xlim([1 ScaleNatural(end)])




% loglog(1:NumSteps, ErrorCutoff, 'k', 'Linewidth', 4)
% hold on
% loglog(1:NumSteps, ErrorQuad1, 'b--', 'Linewidth', 4)
% loglog(1:NumSteps, ErrorQuad2, 'r-.', 'Linewidth', 4)
% xlabel('Number of steps $n$', 'Interpreter','Latex', 'Fontsize', 50)
% xlabel('Number of steps \it n', 'Interpreter', 'Latex', 'Fontsize', 50)


% % % 1st Order Line
% logx = [3,20];
% logy = (30e-3).*logx.^(-1);
% loglog(logx, logy,'k-','linewidth',3)
% 
% text1 = text(4,4*10^(-3),'$1^{st}$ Order','Interpreter','latex');
% set(text1,'Rotation',-20)
% set(text1,'FontSize',50)



% % % 2nd Order Line
% logx = [3,20];
% logy = (30e-3).*logx.^(-2);
% loglog(logx, logy,'k-','linewidth',3)






% MaxFreq = 5;
% ws = ( 0 : StepFreq : MaxFreq)';
% ws2 = ws.^2;
% 
% 
% figure
% j = 10;
% plot(ws,log(H(:,j)),'k.')
% hold on
% plot(ws, log(exp(-0.5*ws2 / ScaleNaturalFreq(j,1)^2)),'b-')
% plot(ws, log(exp(-0.5*ws2 * (2*j*tau))), 'g--')
% 
% plot([1/sqrt(tau), 1/sqrt(tau)], [-70, 0], 'r-')
% xlim([0 MaxFreq])
% ylim([-70 0])
% legend('1st Order Implicit Euler TF','Cutoff Frequency Gaussian', 'Natural Gaussian')
% ylabel('Intensity')
% xlabel('Spatial Frequency')
% 
% 
% 
% log(db3)
% log(exp(-0.5*ws2 / ScaleNaturalFreq(j,1)^2));
% 
% ScaleNaturalFreq(j,1)
% ws(1351)/cut
% 
% open ans
% 
% figure
% plot(ws(1:end-2),diff(diff(log(H(:,40)))),'k.')
% hold on
% plot(ws(1:end-2), zeros(size(ws(1:end-2))))
% plot([1/sqrt(tau), 1/sqrt(tau)], [min(diff(diff(log(H(:,40))))), max(diff(diff(log(H(:,40)))))],'b-')
% 
% [a,b]=min(diff(log(H(:,50))))
% ws(b)













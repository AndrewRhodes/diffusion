% Andrew Rhodes
% WVU
% Dec. 2018
% Compare the scale estimation technique by Fadaifard to our technique
% For journal paper. Fig. 1

close all
clear
clc

NumSteps = 100;
tau = 1;
alpha = 1;

MaxFreq = 2;
StepFreq = 0.0001;


%% Natural Scale Growth Method
ScaleNaturalFreq = zeros(NumSteps,1);
ScaleNatural = zeros(NumSteps,1);

for i = 1 : NumSteps - 1   
    ScaleNatural(i+1,1) = sqrt(2*i*alpha*tau);
    ScaleNaturalFreq(i+1,1) = 1/ScaleNatural(i+1,1);
end




%% Cutoff Frequency Method


CutoffFrequency = zeros(NumSteps,1);
ScaleCutoffFreq = zeros(NumSteps,1);
ScaleCutoff = zeros(NumSteps,1);

db3 = 1/sqrt(2);
cut = sqrt(log(2));
ws = 0 : StepFreq : MaxFreq;
ws2 = (ws.^2)';
NumSample = length(ws);

H = zeros(NumSample, NumSteps);
H(:,1) = ones(NumSample, 1);
h = 1 ./ ( ones(NumSample,1) + tau * ws2);


for i = 1 : NumSteps - 1
    
    H(:,i+1) = H(:,i) .* h;
    % Find the frequency at the cutoff values
    [uH, indH] = unique(H(:,i+1));
    CutoffFrequency(i+1,1) = interp1(uH, ws(indH), db3);
    
    ScaleCutoffFreq(i+1,1) = CutoffFrequency(i+1,1) / cut;   
    ScaleCutoff(i+1,1) =  1 / ScaleCutoffFreq(i+1,1);    
    
end




%% Quadratic Fitting Method (Fadaifard)


ScaleQuadFreq = zeros(NumSteps,1);
ScaleQuad1 = zeros(NumSteps,1);
ScaleQuad2 = zeros(NumSteps,1);

ws = (-2 : StepFreq : 2)';
ws2 = ws.^2;
ws4 = ws2.^2;
sumws4 = sum(ws4);
delta = 1.1;

for i = 1 : NumSteps - 1 
        
    lambda = tau * delta^i;
    
    ScaleQuadFreq(i+1,1) = sqrt( 0.5*sumws4 ...
        ./ sum( (ws2').*sum( log(1 + lambda * ones(i,1) * (ws2')), 1) ,2) );
    
    ScaleQuad1(i+1,1) = 1 / ScaleQuadFreq(i+1,1);
end


ws = (-4 : StepFreq : 4)';
ws2 = ws.^2;
ws4 = ws2.^2;
sumws4 = sum(ws4);

for i = 1 : NumSteps - 1 
        
    lambda = tau * delta^i;
    
    ScaleQuadFreq(i+1,1) = sqrt( 0.5*sumws4 ...
        ./ sum( (ws2').*sum( log(1 + lambda * ones(i,1) * (ws2')), 1) ,2) );
    
    ScaleQuad2(i+1,1) = 1 / ScaleQuadFreq(i+1,1);
end


%% Errors compared to natural scale growth



ErrorCutoff = zeros(NumSteps,1);
ErrorQuad1 = zeros(NumSteps,1);
ErrorQuad2 = zeros(NumSteps,1);


for i = 1 : NumSteps
    
   ErrorCutoff(i,1) =  abs(ScaleNatural(i,1) - ScaleCutoff(i,1)) / ScaleNatural(i,1);
   ErrorQuad1(i,1) =  abs(ScaleNatural(i,1) - ScaleQuad1(i,1)) / ScaleNatural(i,1);
   ErrorQuad2(i,1) =  abs(ScaleNatural(i,1) - ScaleQuad2(i,1)) / ScaleNatural(i,1);
end


% figure
% plot(1:NumSteps, ErrorCutoff, 'b')
% hold on
% plot(1:NumSteps, ErrorQuad, 'r')


% figure
% semilogy(1:NumSteps, ErrorCutoff, 'b')
% hold on
% semilogy(1:NumSteps, ErrorQuad, 'r')



figure
loglog(1:NumSteps, ErrorCutoff, 'k', 'Linewidth', 4)
hold on
loglog(1:NumSteps, ErrorQuad1, 'b--', 'Linewidth', 4)
loglog(1:NumSteps, ErrorQuad2, 'r-.', 'Linewidth', 4)
% xlabel('Number of steps $n$', 'Interpreter','latex', 'Fontsize', 50)
xlabel('Number of steps \it n', 'Interpreter', 'tex', 'Fontsize', 50)
ylabel('Relative Error', 'Fontsize', 50)
ax = gca;
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;

% % % 1st Order Line
logx = [3,20];
logy = (30e-3).*logx.^(-1);
loglog(logx, logy,'k-','linewidth',3)

text1 = text(4,4*10^(-3),'$1^{st}$ Order','Interpreter','latex');
set(text1,'Rotation',-20)
set(text1,'FontSize',50)













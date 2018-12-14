% Andrew Rhodes
% WVU
% September, 2018
% 
% (1) Build Histograms for each individual feature: alpha, phi, beta, theta.
% (2) Match histograms against themselves with nchoosek combinations: alpha,
% alpha-theta, alpha-phi-beta, alpha-phi-beta-theta etc.
% (3) Check the accuracy, repeatability, distance between descriptors. 
%
% t_sigma: ratio between scales used to eliminate features descriptors of
% different scale keypoints.
% 
% t_dist: distance limit between descriptors to say they have matched




% Features = buildHistogram(PointCloud, NMSKeypoint);

% save('Features.mat','Features','-v7.3')
load('Features.mat','Features')

MatchDist = @(A, B) (1 - (1 + sum( min(A, B))) / (1 + sum( max(A, B))) ); % Anonymous Function.


c = 0;
Dist.Data = zeros(Features.Count^2, 18);
WaitBar = waitbar(0, sprintf('Matching Histogram %i of %i', 0, Features.Count));

Dist.Text = {'Feature 1', 'Feature 2', 'Alpha Hist Match', 'Phi Hist Match',...
    'Beta Hist Match', 'Theta Hist Match', 'Alpha,Phi Hist Match',...
    'Alpha,Beta Hist Match', 'Alha,Theta Hist Match', 'Phi,Beta Hist Match',...
    'Phi,Theta Hist Match', 'Beta,Theta Hist Match', 'Alpha,Phi,Beta Hist Match',...
    'Alpha,Phi,Theta Hist Match', 'Alpha,Beta,Theta Hist Match',...
    'Phi,Beta,Theta Hist Match', 'Alpha,Phi,Beta,Theta Hist Match'};

% This loop takes forever and takes 15GB RAM. 
for i = 1 : Features.Count
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform by cellfun
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%     ScaleRatios = Features.Scale / Features.Scale(i);
%     ScaleRatios(ScaleRatios > 1) = 1 ./ (ScaleRatios(ScaleRatios > 1));
%     
%     
%     Dist.Data((i-1)*Features.Count+1:i*Features.Count, 1:18) = [i*ones(Features.Count,1), (1:Features.Count)', ScaleRatios, ...
%         cellfun(MatchDist, mat2cell(repmat(Features.HistAlpha(i,:), Features.Count, 1), ones(Features.Count,1), 45), mat2cell( Features.HistAlpha, ones(Features.Count,1), 45)), ...
%         cellfun(MatchDist, mat2cell(repmat(Features.HistPhi(i,:), Features.Count, 1), ones(Features.Count,1), 45), mat2cell( Features.HistPhi, ones(Features.Count,1), 45)), ...
%         cellfun(MatchDist, mat2cell(repmat(Features.HistBeta(i,:), Features.Count, 1), ones(Features.Count,1), 45), mat2cell( Features.HistBeta, ones(Features.Count,1), 45)), ...
%         cellfun(MatchDist, mat2cell(repmat(Features.HistTheta(i,:), Features.Count, 1), ones(Features.Count,1), 45), mat2cell( Features.HistTheta, ones(Features.Count,1), 45)), ...
%         ...
%         cellfun(MatchDist, mat2cell(repmat([Features.HistAlpha(i,:), Features.HistPhi(i,:)], Features.Count, 1), ones(Features.Count,1), 2*45),...
%         mat2cell( [Features.HistAlpha, Features.HistPhi], ones(Features.Count,1), 2*45)),...
%         ...
%         cellfun(MatchDist, mat2cell(repmat([Features.HistAlpha(i,:), Features.HistBeta(i,:)], Features.Count, 1), ones(Features.Count,1), 2*45),...
%         mat2cell( [Features.HistAlpha, Features.HistBeta], ones(Features.Count,1), 2*45)),...
%         ...
%         cellfun(MatchDist, mat2cell(repmat([Features.HistAlpha(i,:), Features.HistTheta(i,:)], Features.Count, 1), ones(Features.Count,1), 2*45),...
%         mat2cell( [Features.HistAlpha, Features.HistTheta], ones(Features.Count,1), 2*45)),...
%         ...
%         cellfun(MatchDist, mat2cell(repmat([Features.HistPhi(i,:), Features.HistBeta(i,:)], Features.Count, 1), ones(Features.Count,1), 2*45),...
%         mat2cell( [Features.HistPhi, Features.HistBeta], ones(Features.Count,1), 2*45)),...
%         ...
%         cellfun(MatchDist, mat2cell(repmat([Features.HistPhi(i,:), Features.HistTheta(i,:)], Features.Count, 1), ones(Features.Count,1), 2*45),...
%         mat2cell( [Features.HistPhi, Features.HistTheta], ones(Features.Count,1), 2*45)),...
%         ...
%         cellfun(MatchDist, mat2cell(repmat([Features.HistBeta(i,:), Features.HistTheta(i,:)], Features.Count, 1), ones(Features.Count,1), 2*45),...
%         mat2cell( [Features.HistBeta, Features.HistTheta], ones(Features.Count,1), 2*45)),...
%         ...
%         cellfun(MatchDist, mat2cell(repmat([Features.HistAlpha(i,:), Features.HistPhi(i,:), Features.HistBeta(i,:)], Features.Count, 1), ones(Features.Count,1), 3*45),...
%         mat2cell( [Features.HistAlpha, Features.HistPhi, Features.HistBeta], ones(Features.Count,1), 3*45)),...
%         ...
%         cellfun(MatchDist, mat2cell(repmat([Features.HistAlpha(i,:), Features.HistPhi(i,:), Features.HistTheta(i,:)], Features.Count, 1), ones(Features.Count,1), 3*45),...
%         mat2cell( [Features.HistAlpha, Features.HistPhi, Features.HistTheta], ones(Features.Count,1), 3*45)),...
%         ...
%         cellfun(MatchDist, mat2cell(repmat([Features.HistAlpha(i,:), Features.HistBeta(i,:), Features.HistTheta(i,:)], Features.Count, 1), ones(Features.Count,1), 3*45),...
%         mat2cell( [Features.HistAlpha, Features.HistBeta, Features.HistTheta], ones(Features.Count,1), 3*45)),...
%         ...
%         cellfun(MatchDist, mat2cell(repmat([Features.HistPhi(i,:), Features.HistBeta(i,:), Features.HistTheta(i,:)], Features.Count, 1), ones(Features.Count,1), 3*45),...
%         mat2cell( [Features.HistPhi, Features.HistBeta, Features.HistTheta], ones(Features.Count,1), 3*45)),...
%         ...
%         cellfun(MatchDist, mat2cell(repmat([Features.HistAlpha(i,:), Features.HistPhi(i,:), Features.HistBeta(i,:), Features.HistTheta(i,:)], Features.Count, 1), ones(Features.Count,1), 4*45),...
%         mat2cell( [Features.HistAlpha, Features.HistPhi, Features.HistBeta, Features.HistTheta], ones(Features.Count,1), 4*45))];



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform by for-loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%     ScaleRatios = Features.Scale / Features.Scale(i);
%     ScaleRatios(ScaleRatios > 1) = 1 ./ (ScaleRatios(ScaleRatios > 1));
    

%     for j = 1 : Features.Count
%            c = c+1;
%         Dist.Data(c, 1:18) = [i, j, ScaleRatios(j),...
%             ...
%             MatchDist(Features.HistAlpha(i,:), Features.HistAlpha(j,:)), ...
%             MatchDist(Features.HistPhi(i,:), Features.HistPhi(j,:)), ...
%             MatchDist(Features.HistBeta(i,:), Features.HistBeta(j,:)), ... 
%             MatchDist(Features.HistTheta(i,:), Features.HistTheta(j,:)), ...
%             ...
%             MatchDist([Features.HistAlpha(i,:), Features.HistPhi(i,:)], ...
%             [Features.HistAlpha(j,:), Features.HistPhi(j,:)]), ...
%             ...
%             MatchDist([Features.HistAlpha(i,:), Features.HistBeta(i,:)], ...
%             [Features.HistAlpha(j,:), Features.HistBeta(j,:)]), ...
%             ...
%             MatchDist([Features.HistAlpha(i,:), Features.HistTheta(i,:)], ...
%             [Features.HistAlpha(j,:), Features.HistTheta(j,:)]), ...
%             ...
%             MatchDist([Features.HistPhi(i,:), Features.HistBeta(i,:)], ...
%             [Features.HistPhi(j,:), Features.HistBeta(j,:)]), ...
%             ...
%             MatchDist([Features.HistPhi(i,:), Features.HistTheta(i,:)], ...
%             [Features.HistPhi(j,:), Features.HistTheta(j,:)]), ...
%             ...
%             MatchDist([Features.HistBeta(i,:), Features.HistTheta(i,:)], ...
%             [Features.HistBeta(j,:), Features.HistTheta(j,:)]), ...
%             ...            
%             MatchDist([Features.HistAlpha(i,:), Features.HistPhi(i,:), Features.HistBeta(i,:)], ...
%             [Features.HistAlpha(j,:), Features.HistPhi(j,:), Features.HistBeta(j,:)]), ...
%             ...
%             MatchDist([Features.HistAlpha(i,:), Features.HistPhi(i,:), Features.HistTheta(i,:)], ...
%             [Features.HistAlpha(j,:), Features.HistPhi(j,:), Features.HistTheta(j,:)]), ...
%             ...
%             MatchDist([Features.HistAlpha(i,:), Features.HistBeta(i,:), Features.HistTheta(i,:)], ...
%             [Features.HistAlpha(j,:), Features.HistBeta(j,:), Features.HistTheta(j,:)]), ...
%             ...
%             MatchDist([Features.HistPhi(i,:), Features.HistBeta(i,:), Features.HistTheta(i,:)], ...
%             [Features.HistPhi(j,:), Features.HistBeta(j,:), Features.HistTheta(j,:)]), ...
%             ...
%             MatchDist([Features.HistAlpha(i,:), Features.HistPhi(i,:), Features.HistBeta(i,:), Features.HistTheta(i,:)], ...
%             [Features.HistAlpha(j,:), Features.HistPhi(j,:), Features.HistBeta(j,:), Features.HistTheta(j,:)])];
%         
%     end
    
    waitbar(i/Features.Count, WaitBar, sprintf('Matching Histogram %i of %i', i, Features.Count));
end

waitbar(i/Features.Count, WaitBar, sprintf('Matching Histograms Complete'));
close(WaitBar)

% save('Dist.mat','Dist','-v7.3')
load('Dist.mat','Dist')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



t_sigma = 0.7; % 
t_dist = 0.1;


DistMatchLessThanThresh = zeros(Features.Count,15);
DistMatchLessThanThreshScaleLessThenThres = zeros(Features.Count,15);
DistMatchNextClosestScaleThresh = zeros(Features.Count,15);
DistMatchNextClosest = zeros(Features.Count,15);
ClosestDescriptorMatchLocationScaleThresh = zeros(Features.Count,15);
ClosestDescriptorMatchLocation = zeros(Features.Count,15);
DistMatchAllFeaturesLessThanThresh = zeros(Features.Count,5);


WaitBar = waitbar(0, sprintf('Matching Histogram %i of %i', 0, Features.Count));

for i = 1 : Features.Count
    
    ScaleRatios = Features.Scale / Features.Scale(i);
    ScaleRatios(ScaleRatios > 1) = 1 ./ (ScaleRatios(ScaleRatios > 1));
    
	ScalePositions = find(ScaleRatios > t_sigma);
    NumScaleMatches = length(ScalePositions);
    
    CurrentDist = Dist.Data((i-1)*Features.Count + 1 : i*Features.Count,4:18);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Number of Descriptor Match Distances Less Than t_dist
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    DistMatchLessThanThresh(i,:) = sum(CurrentDist < t_dist,1);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Number of Descriptor Match Distances Less Than t_dist
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    DistMatchLessThanThreshScaleLessThenThres(i,:) = sum(CurrentDist(ScalePositions,:) < t_dist,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Number of Descriptor Match Distances Less Than many distances
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DistMatchAllFeaturesLessThanThresh(i,1:5) = [nnz(CurrentDist(:,end)<0.1), nnz(CurrentDist(:,end)<0.2), nnz(CurrentDist(:,end)<0.3), nnz(CurrentDist(:,end)<0.4), nnz(CurrentDist(:,end)<0.5)];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Distance to Next Closest Match using Scale Threshold t_sigma
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [SortVal, SortPos] = sort(CurrentDist(ScalePositions,:), 'ascend');    
    
    DistMatchNextClosestScaleThresh(i,:) = diff(SortVal(1:2,:),[],1);
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Distance to Next Closest Match
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [SortVal, SortPos] = sort(CurrentDist, 'ascend');    
    
    DistMatchNextClosest(i,:) = diff(SortVal(1:2,:),[],1);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Closest Descriptor Matching Location using Scale Thesh t_sigma, should be i
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [MinVal, MinPos] = min(CurrentDist(ScalePositions,:),[],1);
    
    ClosestDescriptorMatchLocationScaleThresh(i,:) = ScalePositions(MinPos)'; % Should be i
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Closest Descriptor Matching Location, should be i
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    [MinVal, MinPos] = min(CurrentDist,[],1);
    ClosestDescriptorMatchLocation(i,:) = MinPos; % Should be i
    
    
    
    waitbar(i/Features.Count, WaitBar, sprintf('Matching Histogram %i of %i', i, Features.Count));
    
    
end

waitbar(i/Features.Count, WaitBar, sprintf('Matching Histograms Complete'));
close(WaitBar)


modelData = {'$\alpha$', '$\phi$', '$\beta$', '$\theta$', '$\alpha\phi$', '$\alpha\beta$',...
    '$\alpha\theta$', '$\phi\beta$', '$\phi\theta$', '$\beta\theta$',...
    '$\alpha\phi\beta$', '$\alpha\phi\theta$', '$\alpha\beta\theta$', '$\phi\beta\theta$', '$\alpha\phi\beta\theta$'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Number of Descriptor Match Distances Less Than t_dist


figure
hold on
Data = [];
Points = [];
for j = 1 : 15
    Data = [Data; DistMatchLessThanThresh(:,j)];
    Points = [Points; j*ones(Features.Count,1)];
end
bp = boxplot(Data, Points);
set(bp,'Linewidth',2,'Markersize',8)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XTickLabel = modelData;
ax.YAxis.FontSize = 50;
ax.YScale = 'log';
ax.XAxis.FontSize = 30;
xtickangle(ax,-30)
xlabel('Feature Combinations','Fontsize',50)
ylabel('Number','Fontsize',50)
title(sprintf('Descriptor Matches, $t_{dist} < %0.1f$',t_dist), 'Fontsize',50, 'interpreter','latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Number of Descriptor Match Distances Less Than t_dist, and scale less
%     than t_scale


figure
hold on
Data = [];
Points = [];
for j = 1 : 15
    Data = [Data; DistMatchLessThanThreshScaleLessThenThres(:,j)];
    Points = [Points; j*ones(Features.Count,1)];
end
bp = boxplot(Data, Points);
set(bp,'Linewidth',2,'Markersize',8)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XTickLabel = modelData;
ax.YAxis.FontSize = 50;
ax.YScale = 'log';
ax.XAxis.FontSize = 30;
xtickangle(ax,-30)
xlabel('Feature Combinations','Fontsize',50)
ylabel('Number','Fontsize',50)
title(sprintf('Descriptor Matches, $t_{dist} < %0.1f$, $t_{scale} > %0.1f$', t_dist, t_sigma), 'Fontsize',50, 'interpreter','latex')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Number of Descriptor Match Distances Less Than many distances

figure
hold on
Data = [];
Points = [];
for j = 1 : 5
    Data = [Data; DistMatchAllFeaturesLessThanThresh(:,j)];
    Points = [Points; j*ones(Features.Count,1)];
end
bp = boxplot(Data, Points);
set(bp,'Linewidth',2,'Markersize',8)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XTickLabel = {'0.1','0.2','0.3','0.4','0.5'};
ax.YAxis.FontSize = 50;
ax.XAxis.FontSize = 50;
xlabel(sprintf('Values of $t_\\sigma$'),'FontSize', 50, 'interpreter','latex')
ylabel('Number','FontSize', 50)
title(sprintf('$\\alpha\\phi\\beta\\theta$ Descriptor Matches'),'FontSize', 50, 'interpreter','latex')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Distance to Next Closest Match using Scale Threshold t_sigma


figure
hold on
Data = [];
Points = [];
for j = 1 : 15
    Data = [Data; DistMatchNextClosestScaleThresh(:,j)];
    Points = [Points; j*ones(Features.Count,1)];
end
bp = boxplot(Data, Points);
set(bp,'Linewidth',2,'Markersize',8)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XTickLabel = modelData;
ax.YAxis.FontSize = 50;
ax.XAxis.FontSize = 30;
xtickangle(ax,-30)
xlabel('Feature Combinations','Fontsize',50)
ylabel('Match Distance','Fontsize',50)
title(sprintf('Next Closest Distance, $t_{\\sigma}$ = %.1f', t_sigma),'Fontsize',50, 'interpreter', 'latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Distance to Next Closest Match

figure
hold on
Data = [];
Points = [];
for j = 1 : 15
    Data = [Data; DistMatchNextClosest(:,j)];
    Points = [Points; j*ones(Features.Count,1)];
end
bp = boxplot(Data, Points);
set(bp,'Linewidth',2,'Markersize',8)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XTickLabel = modelData;
ax.YAxis.FontSize = 50;
ax.XAxis.FontSize = 30;
xtickangle(ax,-30)
xlabel('Feature Combinations','Fontsize',50)
ylabel('Match Distance','Fontsize',50)
title(sprintf('Next Closest Distance'),'Fontsize',50)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Closest Descriptor Matching Location using Scale Thesh t_sigma, should be i


MatchCorrect = ClosestDescriptorMatchLocationScaleThresh - repmat((1:Features.Count)',1,15);
Error = zeros(15,1);
for i = 1 : 15
    Error(i) = (Features.Count - nnz(MatchCorrect(:,i))) / Features.Count;
end


figure
bar(Error)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XTickLabel = modelData;
ax.YAxis.FontSize = 50;
ax.XAxis.FontSize = 30;
xtickangle(ax,-30)
xlabel('Feature Combinations','Fontsize',50)
ylabel('Percentage','Fontsize',50)
title(sprintf('Repeat Selection of Self, $t_{\\sigma}$ = %.1f', t_sigma), 'Fontsize',50, 'interpreter', 'latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Closest Descriptor Matching Location, should be i


MatchCorrect = ClosestDescriptorMatchLocation - repmat((1:Features.Count)',1,15);
Error = zeros(15,1);
for i = 1 : 15
    Error(i) = (Features.Count - nnz(MatchCorrect(:,i))) / Features.Count;
end


figure
bar(Error)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XTickLabel = modelData;
ax.YAxis.FontSize = 50;
ax.XAxis.FontSize = 30;
xtickangle(ax,-30)
xlabel('Feature Combinations','Fontsize',50)
ylabel('Percentage','Fontsize',50)
title('Repeat Selection of Self', 'Fontsize',50)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Example Feature Descriptor using alpha-phi-beta-theta

[a,b]=max(Features.Scale);

figure
plot([Features.HistAlpha(b,:), Features.HistPhi(b,:), Features.HistBeta(b,:), Features.HistTheta(b,:)],'linewidth',4)
xlabel(sprintf('Features $\\alpha, \\phi, \\beta, \\theta$. 45 bins each'),'interpreter', 'latex')
ylabel('Percentage')
title('Example Feature Descriptor')
ax = gca;
xlim([1 180])
ax.XTick = [1, 45, 90, 135, 180];
ax.YAxis.FontSize = 50;
ax.XAxis.FontSize = 50;










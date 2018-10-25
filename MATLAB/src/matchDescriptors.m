% Andrew Rhodes
% WVU
% September 2018
%
% Match feature descriptors between model and scene.
%
% input[]:
% input[]:
% input[]:
% input[]:
% 
% output[]:



% % % % Correspondences.Features = [Scene, Model]; % % % %



function Correspondences = matchDescriptors(SceneFeatures, ModelFeatures, t_scale, t_dist)


% Anonymous Function to match Feature Descriptors
DescriptorMatchFunction = @(A, B) (1 - (1 + sum( min(A, B))) / (1 + sum( max(A, B))) ); 


% t_scale = 0.7;
% t_dist = 0.2;
sp = 0;
ep = 0;

% % % % Correspondences.Features = [Scene, Model]; % % % %
Correspondences.Features = zeros(SceneFeatures.Count*ModelFeatures.Count, 2);
Correspondences.Distance = zeros(SceneFeatures.Count*ModelFeatures.Count, 1);


WaitBar = waitbar(0, sprintf('Matching Histogram %i of %i', 0, SceneFeatures.Count));

for i = 1 : SceneFeatures.Count
    
    % Limit the matches by scale similarity
    ScaleRatios = min(ModelFeatures.Scale, SceneFeatures.Scale(i)) ./ max(ModelFeatures.Scale, SceneFeatures.Scale(i));
    
    ScalePositions = find(ScaleRatios > t_scale);  
    NumScaleMatches = length(ScalePositions);
    
    
    
    sp = ep + 1;
    ep = sp + NumScaleMatches - 1;
    
    Correspondences.Features(sp:ep,1:2) = [i*ones(NumScaleMatches,1), ScalePositions];
    
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %     Match Descriptors
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%         
%     %     Match using cellfun
% %     MatchDistance = cellfun(DescriptorMatchFunction, mat2cell(repmat([SceneFeatures.HistAlpha(i,:), SceneFeatures.HistPhi(i,:), SceneFeatures.HistBeta(i,:), SceneFeatures.HistTheta(i,:)], NumScaleMatches, 1), ones(NumScaleMatches,1), 4*45),...
% %         mat2cell( [ModelFeatures.HistAlpha(ScalePositions,:), ModelFeatures.HistPhi(ScalePositions,:), ModelFeatures.HistBeta(ScalePositions,:), ModelFeatures.HistTheta(ScalePositions,:)], ones(NumScaleMatches,1), 4*45));
% 
% %     Match using for-loop
%     MatchDistance = zeros(NumScaleMatches,1);
%     for j = 1 : NumScaleMatches
%         
%         MatchDistance(j,1) = DescriptorMatchFunction([SceneFeatures.HistAlpha(i,:), SceneFeatures.HistPhi(i,:), SceneFeatures.HistBeta(i,:), SceneFeatures.HistTheta(i,:)],...
%             [ModelFeatures.HistAlpha(ScalePositions(j),:), ModelFeatures.HistPhi(ScalePositions(j),:), ModelFeatures.HistBeta(ScalePositions(j),:), ModelFeatures.HistTheta(ScalePositions(j),:)]);
%         
%     end
%     
%     % Limit matches to less than matching distance threshold
%     MatchDistanceThreshold = MatchDistance < t_dist;
%     NumMatchDistanceThreshold = nnz(MatchDistanceThreshold);
%     
%     if NumMatchDistanceThreshold
%         % This scene descriptor matches a model descriptor!!
%         sp = ep + 1;
%         ep = sp + NumMatchDistanceThreshold - 1;
%         
%         Correspondences.Features(sp:ep,1:2) = [i*ones(NumMatchDistanceThreshold,1), ScalePositions(MatchDistanceThreshold)];
% %         Correspondences.Features(sp:ep,1:2) = [SceneFeatures.LocationIndex(i,1)*ones(NumMatchDistanceThreshold,1), ModelFeatures.LocationIndex(ScalePositions(MatchDistanceThreshold))];
%         
%         Correspondences.Distance(sp:ep,1) = MatchDistance(MatchDistanceThreshold);
%     end
    
    
    
    waitbar(i/SceneFeatures.Count, WaitBar, sprintf('Matching Histogram %i of %i', i, SceneFeatures.Count));
end


waitbar(i/SceneFeatures.Count, WaitBar, sprintf('Matching Histograms Complete'));

close(WaitBar)


ZeroLogic = (Correspondences.Features(:,1) == 0) & (Correspondences.Features(:,2) == 0);
Correspondences.Features(ZeroLogic,:) = [];
Correspondences.Distance(ZeroLogic,:) = [];
Correspondences.Count = length(Correspondences.Distance);

end






% function Correspondences = matchDescriptors(ModelFeatures, SceneFeatures, t_scale, t_dist)
% 
% 
% % Anonymous Function to match Feature Descriptors
% DescriptorMatchFunction = @(A, B) (1 - (1 + sum( min(A, B))) / (1 + sum( max(A, B))) ); 
% 
% 
% % t_scale = 0.7;
% % t_dist = 0.2;
% sp = 0;
% ep = 0;
% 
% % % % % Correspondences.Features = [Scene, Model]; % % % %
% Correspondences.Features = zeros(SceneFeatures.Count*ModelFeatures.Count, 2);
% Correspondences.Distance = zeros(SceneFeatures.Count*ModelFeatures.Count, 1);
% 
% 
% WaitBar = waitbar(0, sprintf('Matching Histogram %i of %i', 0, ModelFeatures.Count));
% 
% 
% for i = 1 : ModelFeatures.Count
%     
%     ScaleRatio = SceneFeatures.Scale ./ ModelFeatures.Scale(i);
%     ScaleRatio( ScaleRatio > 1 ) = 1 ./ ScaleRatio( ScaleRatio > 1 );
%     
%     ScaleCriterion = find(ScaleRatio > t_scale);
%     NumScaleMatches = length(ScaleCriterion);
%     
%     MatchDist = zeros(NumScaleMatches,1);
%     
%     for j = 1 : NumScaleMatches
%         
%         ModelHist = [ModelFeatures.HistAlpha(i,:), ModelFeatures.HistPhi(i,:), ModelFeatures.HistBeta(i,:), ModelFeatures.HistTheta(i,:)];
%         SceneHist = [SceneFeatures.HistAlpha(ScaleCriterion(j),:), SceneFeatures.HistPhi(ScaleCriterion(j),:), SceneFeatures.HistBeta(ScaleCriterion(j),:), SceneFeatures.HistTheta(ScaleCriterion(j),:)];
%         
%         MatchDist(j,1) = DescriptorMatchFunction(SceneHist, ModelHist);
%         
%     end
%     
%     DistCriterion = (MatchDist < t_dist);
%     
%     NumDistMatches = nnz(DistCriterion);
%     
%     if NumDistMatches
%         sp = ep + 1;
%         ep = sp + NumDistMatches - 1;
%         Correspondences.Features(sp:ep, 1:2) = [ ScaleCriterion(DistCriterion), i*ones(NumDistMatches,1)];
%         Correspondences.Distance(sp:ep, 1) = MatchDist(DistCriterion);
%         
%     end
%     
%     
%  waitbar(i/ModelFeatures.Count, WaitBar, sprintf('Matching Histogram %i of %i', i, ModelFeatures.Count));
% end
% 
% 
% waitbar(i/ModelFeatures.Count, WaitBar, sprintf('Matching Histograms Complete'));
% 
% close(WaitBar)
% 
% 
% ZeroLogic = (Correspondences.Features(:,1) == 0) & (Correspondences.Features(:,2) == 0);
% Correspondences.Features(ZeroLogic,:) = [];
% Correspondences.Distance(ZeroLogic,:) = [];
% Correspondences.Count = length(Correspondences.Distance);
% 
% 
% 
% 
% 
% end

































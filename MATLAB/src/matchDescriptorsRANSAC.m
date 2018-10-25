% Andrew Rhodes
% October 2018
% 
% Use ransac to match the minimum set (3) unique features and estimate pose




% % % % Correspondences.Features = [Scene, Model]; % % % %


function [GroupPose, BestInlierRatio, BestCorrespond, BestPose] = matchDescriptorsRANSAC(Correspondences, SceneFeatures, ModelFeatures, ScenePointCloud, ModelPointCloud)


% NumIter = nchoosek(Correspondences.Count, 3); % Huge, don't actually use

ModelPointCloutLocationTree = KDTreeSearcher(ModelPointCloud.Location, 'Distance', 'Euclidean');

Number = 3;
z = 0.99;
w = Number / Correspondences.Count;
k = log(1-z) / log(1-w^Number); % Initial guess of number of iterations


CurrentIter = 0;
c = 0;

WaitBar = waitbar(0, sprintf('RANSAC %i of %6.f, Inlier %%: %.4f', 0, k, w));

while CurrentIter < k
   
    CurrentIter = CurrentIter + 1;
    loop = 1;

    while loop
 
        CurrentSelection = datasample(1:Correspondences.Count, 3, 'Replace', false);
        CurrentFeatures = Correspondences.Features(CurrentSelection,:);
        
        if (length(unique(CurrentFeatures(:,1)))~=3) || (length(unique(CurrentFeatures(:,2)))~=3)
            loop = 1;
        else
            loop = 0;
        end
    end
    
    c = c + 1;
    GroupPose{c,1} = findPose(CurrentFeatures, SceneFeatures, ModelFeatures);
    
    
%     Inliers = findInliers(GroupPose{c,1}, Correspondences, SceneFeatures, ModelFeatures, ModelPointCloutLocationTree);
    Inliers = findInliers(GroupPose{c,1}, ScenePointCloud, ModelPointCloutLocationTree);
    
    if Inliers.Count > 0
       w2 = Inliers.Count / Correspondences.Count;
       
       if w2 > w
           w = w2;
           k = 2*log(1-z) / log(1-w^Number);
           BestInlierRatio = w;
           BestCorrespond = Inliers.LocationIndex;
           BestPose = GroupPose{c,1};
       end
       
    end
    
    waitbar(CurrentIter/k, WaitBar, sprintf('RANSAC %i of %6.f, Inlier %%: %.4f', CurrentIter, k, w));
    
    
end


waitbar(1, WaitBar, 'RANSAC Complete');

close(WaitBar)






end




% function Inliers = findInliers(Pose, Correspondences, SceneFeatures, ModelFeatures, ModelFeatureLocationTree)
% R_s2m = Pose.Rotation'; % Scene to Model
% t_s2m = - (Pose.Rotation' * Pose.translation')'; % Scene to Model
% 
% SceneLocation2ModelLocation = bsxfun(@plus, (R_s2m * SceneFeatures.Location')', t_s2m);
% 
% 
% CorrInDist = sqrt(sum((ModelFeatures.Location(Correspondences.Features(:,2),:) ...
%              - SceneLocation2ModelLocation(Correspondences.Features(:,1),:)).^2,2))...
%              < SceneFeatures.Scale(Correspondences.Features(:,1));
% 
% 
% % searcher = @(Tree, Point, Radius) rangesearch(Tree, Point, Radius);
% % 
% % Position = cellfun(searcher, repmat({ModelFeatureLocationTree},SceneFeatures.Count,1), mat2cell(SceneFeatures.Location, ones(SceneFeatures.Count,1), 3), num2cell(3*SceneFeatures.Scale));
% 
% % if any(~cellfun(@isempty, Position))
% if nnz(CorrInDist)  
%    Inliers.LocationIndex = find(CorrInDist);
%    Inliers.Count = nnz(CorrInDist);
% else
%     Inliers.Count = 0;
% end
% 
% end


function Inliers = findInliers(Pose, ScenePointCloud, ModelPointCloutLocationTree)

R_s2m = Pose.Rotation'; % Scene to Model
t_s2m = - (Pose.Rotation' * Pose.translation')'; % Scene to Model

SceneLocation2ModelLocation = bsxfun(@plus, (R_s2m * ScenePointCloud.Location')', t_s2m);

% searcher = @(Tree, Point, Radius) rangesearch(Tree, Point, Radius);
% 
% Position = cellfun(searcher, repmat({ModelPointCloutLocationTree},ScenePointCloud.LocationCount,1),...
%     mat2cell(SceneLocation2ModelLocation, ones(ScenePointCloud.LocationCount,1), 3), num2cell(2*ScenePointCloud.Resolution*ones(ScenePointCloud.LocationCount,1)));

Position = rangesearch(ModelPointCloutLocationTree, SceneLocation2ModelLocation, 2*ScenePointCloud.Resolution);

CellNotEmtpy = ~cellfun(@isempty, Position);

if any(CellNotEmtpy)
   Inliers.LocationIndex = find(CellNotEmtpy);
   Inliers.Count = nnz(CellNotEmtpy);
else
    Inliers.Count = 0;
end

end









function Pose = findPose(CurrentFeatures, SceneFeatures, ModelFeatures)



y = SceneFeatures.Location(CurrentFeatures(:,1),:);
ybar = sum(y,1) / size(y,1);
ycheck = bsxfun(@minus, y, ybar);

p = ModelFeatures.Location(CurrentFeatures(:,2),:);
pbar = sum(p,1) / size(p,1);
pcheck = bsxfun(@minus, p, pbar);

% B = ycheck' * pcheck;

B = SceneFeatures.Normal(CurrentFeatures(:,1),:)' * ModelFeatures.Normal(CurrentFeatures(:,2),:);

[U, ~, V] = svd(B);

R = U * [1, 0, 0;
    0, 1, 0;
    0, 0, det(U)*det(V)] * V';

R(abs(R)<eps) = 0;
t = ybar - (R * pbar')';
t(abs(t)<eps) = 0;

Pose.Rotation = R;
Pose.translation = t;

        
        
        
        
        
end























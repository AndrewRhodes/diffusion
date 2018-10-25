% Andrew Rhodes
% WVU
% September 2018
%
% find rigid rotation and translation from model to scene
% This is a first estimate and is updated with ICP later





function GroupPose = findRigidPose(GroupCorrespondence, GroupValues, ModelFeatures, SceneFeatures)


NumGroupCorrespondence = length(GroupCorrespondence);

GroupPose = cell(NumGroupCorrespondence,1);
RemoveGroup = [];
WaitBar = waitbar(0, sprintf('Finding Pose %i of %i', 0, NumGroupCorrespondence));

for i = 1 : NumGroupCorrespondence
    
    if size(GroupCorrespondence{i,1},1) >= 3
        
        y = SceneFeatures.Location(GroupCorrespondence{i,1}(:,1),:);
        ybar = sum(y,1) / size(y,1);
        ycheck = bsxfun(@minus, y, ybar);
        
        p = ModelFeatures.Location(GroupCorrespondence{i,1}(:,2),:);
        pbar = sum(p,1) / size(p,1);
        pcheck = bsxfun(@minus, p, pbar);
        
        B = ycheck' * pcheck;
        
        B = SceneFeatures.Normal(GroupCorrespondence{i,1}(:,1),:)' * ModelFeatures.Normal(GroupCorrespondence{i,1}(:,2),:);


        [U, ~, V] = svd(B);
        
        R = U * [1, 0, 0;
            0, 1, 0;
            0, 0, det(U)*det(V)] * V';
        
        R(abs(R)<eps) = 0;
        t = ybar - (R * pbar')';
        t(abs(t)<eps) = 0;
        
        GroupPose{i,1}.Rotation = R;
        GroupPose{i,1}.translation = t;
        GroupPose{i,1}.weight = sum(GroupValues{i,1});
        
    else
        RemoveGroup = [RemoveGroup; i];
    end
    
    waitbar(i/NumGroupCorrespondence, WaitBar, sprintf('Finding Pose %i of %i', i, NumGroupCorrespondence));
    
    
end

waitbar(1, WaitBar, 'Finding Pose Complete');

close(WaitBar)

GroupPose(RemoveGroup) = [];


end






% v = [0,2,1];
% v = v/norm(v);
% e = 25*pi/180;
% Rot = v_to_R(e*v)
%
% trans = [80,50,50]
%
%
% OldLoc = SceneFeatures.Location;
%
% SceneFeatures.Location = bsxfun(@plus, (Rot*OldLoc')', trans);











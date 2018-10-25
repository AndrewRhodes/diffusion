% Andrew Rhodes
% WVU
% September 2018
%
% Group matched descriptors based on geometric consistency.
%
% input[]:
% input[]:
% input[]:
% input[t_gc]: Geometric consistency distance threshold. defaul 0.25
% input[t_corr]: quantity of geometric correspondence threshold. default 0.25
% 
% output[]:


% % % % Correspondences.Features = [Scene, Model]; % % % %




function [GeometricCorrespondence, GroupCorrespondence, GroupValues] = groupKeypoints(Correspondences, SceneFeatures, ModelFeatures, t_gc, t_corr)


% Set in Andrew Johnson (Spin Image) Dissertation. p.41
% t_gc = 0.25;
% t_corr = 0.25;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check distance for Geometric Consistency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set in Andrew Johnson (Spin Image) Dissertation p. 40, Eq. 3.1

c = 0;

WaitBar = waitbar(0, sprintf('Checking Geometric Consistency %i of %i', 0, Correspondences.Count));

for i = 1 : Correspondences.Count
    
    C1 = Correspondences.Features(i,1:2);
    C2 = Correspondences.Features(:,1:2);
    C2(i,:) = [];
    
    
    % C1Model to C2Model
    [CosDistModel, SinDistModel] = findSinCosDist(...
        ModelFeatures.Location(C1(1,2),:),...
        ModelFeatures.Normal(C1(1,2),:),...
        ModelFeatures.Location(C2(:,2),:) );
    
    % C1Scene to C2Scene
    [CosDistScene, SinDistScene] = findSinCosDist(...
        SceneFeatures.Location(C1(1,1),:),...
        SceneFeatures.Normal(C1(1,1),:),...
        SceneFeatures.Location(C2(:,1),:) );
    
    dgc(:,1) = finddgc(CosDistModel, SinDistModel, CosDistScene, SinDistScene);
    

    % C2Model to C1Model
    [CosDistModel, SinDistModel] = findSinCosDist(...
        ModelFeatures.Location(C2(:,2),:),...
        ModelFeatures.Normal(C2(:,2),:),...
        ModelFeatures.Location(C1(1,2),:) );
        
    % C2Scene to C1Scene
    [CosDistScene, SinDistScene] = findSinCosDist(...
        SceneFeatures.Location(C1(:,1),:),...
        SceneFeatures.Normal(C1(:,1),:),...
        SceneFeatures.Location(C2(1,1),:) );
    
    dgc(:,2) = finddgc(CosDistModel, SinDistModel, CosDistScene, SinDistScene);
    
    Dgc = max(dgc, [], 2);
    
    clear dgc
    
    if nnz(Dgc < t_gc) %>= t_corr * Correspondences.Count
        c = c + 1;
        GeometricCorrespondence(c, 1:2) = C1;
        
    end
    
    waitbar(i/Correspondences.Count, WaitBar, sprintf('Checking Geometric Consistency %i of %i', i, Correspondences.Count));
    
    
end

waitbar(1, WaitBar, 'Checking Geometric Consistency Complete');

close(WaitBar)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check Keypoints for Grouping Criterion, Wgc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set in Andrew Johnson (Spin Image) Dissertation p. 41, Eq. 3.2


NumGeometricCorrespondence = size(GeometricCorrespondence, 1);
GroupCorrespondence = cell(NumGeometricCorrespondence,1);
GroupValues = cell(NumGeometricCorrespondence,1);


WaitBar = waitbar(0, sprintf('Grouping Correspondences %i of %i', 0, NumGeometricCorrespondence));

for i = 1 : NumGeometricCorrespondence

    loop = 1;
%     CurrentCorrespondence = zeros(25,2);
    CurrentCorrespondence = GeometricCorrespondence(i,:);
    CurrentValues = []; % zeros(25,1); % 

    
    RemainingCorrespondence = GeometricCorrespondence;
%     RemainingCorrespondence(i,:) = [];
    AllIndex = i;
    
    while loop
        
%         Wgc = zeros(size(RemainingCorrespondence,1), size(CurrentCorrespondence,1));
        wgc = zeros(size(RemainingCorrespondence,1), 2);
        
%         for j = 1 : size(CurrentCorrespondence,1)
%         for j = 1 : c
%             j
            % C1Model to C2Model
            [CosDistModel, SinDistModel] = findSinCosDist(...
                ModelFeatures.Location(CurrentCorrespondence(end,2),:),...
                ModelFeatures.Normal(CurrentCorrespondence(end,2),:),...
                ModelFeatures.Location(RemainingCorrespondence(:,2),:));
            
            % C1Scene to C2Scene
            [CosDistScene, SinDistScene] = findSinCosDist(...
                SceneFeatures.Location(CurrentCorrespondence(end,1),:),...
                SceneFeatures.Normal(CurrentCorrespondence(end,1),:),...
                SceneFeatures.Location(RemainingCorrespondence(:,1),:));
            
            wgc(:,1) = findwgc(CosDistModel, SinDistModel,...
                CosDistScene, SinDistScene,...
                ModelFeatures.Resolution, SceneFeatures.Resolution);
            
   
            
            
            
            % C2Model to C1Model
            [CosDistModel, SinDistModel] = findSinCosDist(...
                ModelFeatures.Location(RemainingCorrespondence(:,2),:),...
                ModelFeatures.Normal(RemainingCorrespondence(:,2),:),...
                ModelFeatures.Location(CurrentCorrespondence(end,2),:));
            
            % C2Scene to C1Scene
            [CosDistScene, SinDistScene] = findSinCosDist(...
                SceneFeatures.Location(RemainingCorrespondence(:,1),:),...
                SceneFeatures.Normal(RemainingCorrespondence(:,1),:),...
                SceneFeatures.Location(CurrentCorrespondence(end,1),:));
            
            wgc(:,2) = findwgc(CosDistModel, SinDistModel,...
                CosDistScene, SinDistScene,...
                ModelFeatures.Resolution, SceneFeatures.Resolution);
                        
            Wgc = max(wgc, [], 2);
            
%             clear wgc
%         end
        
        if ~isempty(AllIndex)
           Wgc(AllIndex,:) = inf; 
        end
        
        [WgcValue, index] = min(max(Wgc,[],2),[],1);
        AllIndex = [AllIndex; index];
        
        %         [WgcValue, index] = min(Wgc,[],1);
%                 WgcValue
        
%         if (c > 25) || (WgcValue > t_gc*t_gc)
%             
%             loop = 0;
%             CurrentValues(CurrentCorrespondence(:,1) == 0, :) = [];
%             CurrentCorrespondence(CurrentCorrespondence(:,1) == 0, :) = [];
%             
%             GroupCorrespondence{i,1} = CurrentCorrespondence;
% %             size(CurrentCorrespondence)
%             
%             GroupValues{i,1} = CurrentValues;
%             
%         else
        if WgcValue < t_gc
            
            loop = 1;
%             c = c + 1;
            CurrentCorrespondence = [CurrentCorrespondence; RemainingCorrespondence(index,:)];
%             CurrentCorrespondence(c,:) = RemainingCorrespondence(index,:);
%             CurrentValues(c,1) = WgcValue;
            CurrentValues = [CurrentValues; WgcValue];
            
%             RemainingCorrespondence(index,:) = [];
            %             size(RemainingCorrespondence)
%             
        else
            loop = 0;
            % %
            %             CurrentCorrespondence(CurrentCorrespondence(:,1) == 0, :) = [];
            %
            [UniqueValues, UniqueIndex, UniqueToCurrentCorrespondences] = unique(CurrentCorrespondence(:,1));
            KeepColumn1 = false(size(CurrentCorrespondence,1),1);
            KeepColumn1(UniqueIndex) = 1;
            
            [UniqueValues, UniqueIndex, UniqueToCurrentCorrespondences] = unique(CurrentCorrespondence(:,2));
            KeepColumn2 = false(size(CurrentCorrespondence,1),1);
            KeepColumn2(UniqueIndex) = 1;
            
            GroupCorrespondence{i,1} = CurrentCorrespondence(KeepColumn1 & KeepColumn2,:);
            
%             GroupCorrespondence{i,1} = CurrentCorrespondence;

            %             size(CurrentCorrespondence)
            %
            GroupValues{i,1} = CurrentValues;
            %

        end
        
        
        clear wgc Wgc
        
    end
    
    clear CurrentCorrespondence CurrentValues AllIndex
    
    waitbar(i/NumGeometricCorrespondence, WaitBar, sprintf('Grouping Correspondences %i of %i', i, NumGeometricCorrespondence));
    
    
end

waitbar(1, WaitBar, 'Grouping Correspondences Complete');

close(WaitBar)



end





function wgc = findwgc(CosDistModel, SinDistModel, CosDistScene, SinDistScene, ModelResolution, SceneResolution)
% Find the group consistency weighting between model and scene
% 
% input[CosDistModel]: 
% input[SinDistModel]:
% input[CosDistScene]:
% input[SinDistScene]:
% input[ModelResolution]: Average edge length of model
% input[SceneResolution]: Average edge length of scene
% 
% output[wgc]: size[n,1]

dgc = finddgc(CosDistModel, SinDistModel, CosDistScene, SinDistScene);

Den = (1 - exp(-sqrt(CosDistModel.^2 + SinDistModel.^2) / (2*4*ModelResolution/2)))...
      .* (1 - exp(-sqrt(CosDistScene.^2 + SinDistScene.^2) / (2*4*SceneResolution/2)));
  
      
% Den = 1 - exp(-sqrt(CosDistModel.^2 + SinDistModel.^2) / (2*4*ModelResolution))...
%        .* exp(-sqrt(CosDistScene.^2 + SinDistScene.^2) / (2*4*SceneResolution))
  
wgc = bsxfun(@rdivide, dgc, Den);

end






function dgc = finddgc(CosDistModel, SinDistModel, CosDistScene, SinDistScene)
% Find the geometric consistency distance between model and scene
% 
% input[CosDistModel]: Cos Distance on Model. size[n,1]
% input[SinDistModel]: Sin Distance on Model. size[n,1]
% input[CosDistScene]: Cos Distance on Scene. size[n,1]
% input[SinDistScene]: Sin Distance on Scene. size[n,1]
% output[dgc]: geometric consistency distance. size[n,1]


dgc = 2 * bsxfun(@rdivide, sqrt(sum(([CosDistModel, SinDistModel] - [CosDistScene, SinDistScene]).^2,2)), ...
    sqrt(sum(([CosDistModel, SinDistModel]).^2,2)) + sqrt(sum(([CosDistScene, SinDistScene]).^2,2)));

end





function [CosDist, SinDist] = findSinCosDist(Location, Normal, OtherLocation)
% Find the Manhattan distances between two oriented points 
% 
% input[Location]: Point location. size[1,3]
% input[Normal]: point normal vector. size [1,3]
% input[OtherLocation]: other points locations. size[n,3]
% output[CosDist]: Cos Distance. size[n,1]
% output[SinDist]: Sin Distance. size[n,1]

Direction = bsxfun(@minus, Location, OtherLocation);
Distance = sqrt(sum(Direction.^2,2));

CosDirAxis = sum(bsxfun(@times, bsxfun(@rdivide, Direction, Distance), Normal),2);
CosDirAxis = max(-1, min(1, CosDirAxis)); % Ensure in range [-1,1]

CosDist = bsxfun(@times, Distance, CosDirAxis);
SinDist = bsxfun(@times, Distance, sqrt(1-CosDirAxis.^2)); % SinDirAxis = sin = sqrt(1-cos^2) = sqrt(1-CosDirAxis.^2)


end





































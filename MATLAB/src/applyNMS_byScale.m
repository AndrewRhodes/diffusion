% Andrew Rhodes
% Non Maximum Suppression of Keypoints
% 
% input[PointCloud]: 
% input[DoG]: Difference of Gaussian structure. size[PointCloud.LocationCount, NumLevels]
% input[Keypoint]: structure with Location, Level, Scale
% input[t_scale]: scalar
% input[t_DoG]: scalar
% input[t_range]: scalar. the range for rangesearch.
%
% output[NMSKeypoint]: Keypoint structure after removing non maxima




function NMSKeypoint = applyNMS_byScale(PointCloud, DoG, Keypoint, t_scale, t_DoG, t_range)


NumKeypoints = length(Keypoint.Location);

KeypointTree = KDTreeSearcher(PointCloud.Location(Keypoint.Location,:), 'Distance', 'euclidean');

Ind = sub2ind(size(DoG), Keypoint.Location, Keypoint.Level);
KpDoG = DoG(Ind);

[KpScaleSort, KpDoGInd] = sort(Keypoint.Scale, 'ascend');

RemoveInd = [];
RemoveIndAll = [];

KpIndRemain = ones(NumKeypoints,1);
Discontinue = 0;

WaitBar = waitbar(0, sprintf('Checking Keypoint %i of %i', 0, NumKeypoints));

for i = 1 : NumKeypoints
    
    CurrentPoint = PointCloud.Location(Keypoint.Location(KpDoGInd(i)),:);
    CurrentScale = KpScaleSort(i);
    CurrentDoG = KpDoG(KpDoGInd(i));
    CurrentInd = KpDoGInd(i);
    
    srange = t_range * PointCloud.ResolutionLocal(Keypoint.Location(KpDoGInd(i)));
    
    if KpIndRemain(KpDoGInd(i))
        
        [Neigh, Dist] = rangesearch(KeypointTree, CurrentPoint, srange);
        Neigh = Neigh{:};
        Dist = Dist{:};
        
        % Remove the repeated point from searching
        Neigh(Neigh == KpDoGInd(i)) = [];
        Dist(Neigh == KpDoGInd(i)) = [];
        
        if ~isempty(Neigh)
            for j = 1 : length(Neigh)
                if KpIndRemain(Neigh(j))
                    
                    NextScale = Keypoint.Scale(Neigh(j));
                    NextDoG = KpDoG(Neigh(j));
                    
                    % if the scales of the two keypoints are similar
                    if ( min(CurrentScale, NextScale) / max(CurrentScale,NextScale) ) > t_scale
                        
                        if abs(NextDoG) > abs(CurrentDoG)
                            
                            RemoveInd = [RemoveInd, Neigh(j)];
                            
                        elseif abs(NextDoG) <= abs(CurrentDoG)
                            
                            RemoveInd = [RemoveInd, KpDoGInd(i)];
                            Discontinue = 1;
                            
                        else
                            warning('NMS needs another criterion.')
                        end
                        
                        
                    end
                    
                    
                end
                
                if Discontinue
                    Discontinue = 0;
                    break;
                end
                
            end
            
            if ~isempty(RemoveInd)
                KpIndRemain(KpDoGInd(RemoveInd)) = 0;
                KpIndRemain(KpDoGInd(i)) = 0;
                RemoveIndAll = [RemoveIndAll, RemoveInd];
                RemoveInd = [];
            end
        end 
    end
    
    
    
    waitbar(i/NumKeypoints, WaitBar, sprintf('Checking Keypoint %i of %i', i, NumKeypoints));
    
end

waitbar(i/NumKeypoints, WaitBar, sprintf('Checking Keypoints Complete'));
close(WaitBar)


NMSKeypoint = Keypoint;

NMSKeypoint.Location(RemoveIndAll) = [];
NMSKeypoint.Scale(RemoveIndAll) = [];
NMSKeypoint.Level(RemoveIndAll) = [];
% NMSKeypoint.LocationCell(RemoveIndAll) = [];



end

























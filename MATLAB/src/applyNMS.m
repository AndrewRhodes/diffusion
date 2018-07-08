% Andrew Rhodes
% Non Maximum Suppression of Keypoints
% 
% input[PointCloud]: 
% input[DoG]: Difference of Gaussian structure. size[PointCloud.LocationCount, NumLevels]
% input[Keypoint]: structure with Location, Level, Scale
% input[t_scale]: scalar
% input[t_DoG]: scalar
%
% output[NMSKeypoint]: Keypoint structure after removing non maxima




function NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_DoG)

% t_scale = 0.9;
% t_DoG = 0.9;


NumKeypoints = length(Keypoint.Location);

KeypointTree = KDTreeSearcher(PointCloud.Location(Keypoint.Location,:), 'Distance', 'euclidean');

Ind = sub2ind(size(DoG), Keypoint.Location, Keypoint.Level);
KpDoG = DoG(Ind);

[KpDoGSort, KpDoGInd] = sort(abs(KpDoG),'descend');

RemoveInd = [];
RemoveIndAll = [];

KpIndRemain = ones(NumKeypoints,1);



for i = 1 : NumKeypoints

    CurrentPoint = PointCloud.Location(Keypoint.Location(KpDoGInd(i)),:);
    CurrentScale = Keypoint.Scale(KpDoGInd(i));
    CurrentDoG = KpDoGSort(i);
    CurrentInd = KpDoGInd(i);
    
    if KpIndRemain(KpDoGInd(i))
        
        [Neigh, Dist] = rangesearch(KeypointTree, CurrentPoint, CurrentScale);
        
        Neigh = Neigh{:};
        Dist = Dist{:};
        
        if length(Neigh) > 1
            for j = 1 : length(Neigh)
                if KpIndRemain(Neigh(j))

                    NextScale = Keypoint.Scale(Neigh(j));
                    NextDoG = KpDoG(Neigh(j));

                    if min(CurrentScale,NextScale)/max(CurrentScale,NextScale) > t_scale
  
                        if abs(NextDoG) > abs(CurrentDoG)
                            % CurrentPoint and NextPoint are in different
                            % locations, with similar scale, but the Dog of
                            % NextPoint is absolutely greater than the DoG of
                            % CurrentPoint. Remove CurrentPoint
                            RemoveInd = [RemoveInd, KpDoGInd(i)];
                                                        
                        elseif abs(NextDoG) < abs(CurrentDoG)
                            % CurrentPoint and NextPoint are in different
                            % locations, with similar scale, but the DoG of
                            % NextPoint is absolutely less than the DoG of
                            % CurrentPoint. Remove NextPoint
                            
                            RemoveInd = [RemoveInd, Neigh(j)];
                            
                        elseif ( (min(abs(NextDoG),abs(CurrentDoG)) / max(abs(NextDoG),abs(CurrentDoG))) > t_DoG ) && (KpDoGInd(i) ~= Neigh(j))
                            % CurrentPoint and NextPoint are in different
                            % locations, with the same DoG and the similar scale.
                            % Randomly remove one of them.
                            
                            RemoveInd = [RemoveInd, randsample([KpDoGInd(i), Neigh(j)],1)];
                            
                        elseif ( (min(abs(NextDoG),abs(CurrentDoG)) / max(abs(NextDoG),abs(CurrentDoG))) > t_DoG ) && (KpDoGInd(i) == Neigh(j))
                            % The Current Point is the Next Point
                            % Nothing Happens
                        else
                            warning('NMS needs another criterion')
                        end
                    end
                end
            end
            
            if ~isempty(RemoveInd)
                KpIndRemain(KpDoGInd(RemoveInd)) = 0;
                KpIndRemain(KpDoGInd(i)) = 0;
                RemoveIndAll = [RemoveIndAll, RemoveInd];
                KpDoGInd(i);
                RemoveInd = [];
            end
        end
    end
    
end


NMSKeypoint = Keypoint;

NMSKeypoint.Location(RemoveIndAll) = [];
NMSKeypoint.Scale(RemoveIndAll) = [];
NMSKeypoint.Level(RemoveIndAll) = [];



end







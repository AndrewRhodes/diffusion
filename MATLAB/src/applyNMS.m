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




function NMSKeypoint = applyNMS(PointCloud, DoG, Keypoint, t_scale, t_range, DoGNormalize, CompareMethod)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove Keypoints that are too small
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RemoveTooSmall = (Keypoint.Scale < PointCloud.ResolutionLocal(Keypoint.LocationIndex));

if isfield(Keypoint, 'Count')
    Keypoint = rmfield(Keypoint, 'Count');
end
FNames = fieldnames(Keypoint);
for i = 1 : length(FNames)
    Keypoint.(FNames{i})(RemoveTooSmall,:) = [];
end
Keypoint.Count = length(Keypoint.LocationIndex);

% Keypoint.Scale(RemoveTooSmall,:) = [];
% Keypoint.Location(RemoveTooSmall,:) = [];
% Keypoint.LocationIndex(RemoveTooSmall,:) = [];
% Keypoint.Level(RemoveTooSmall,:) = [];
% Keypoint.Normal(RemoveTooSmall,:) = [];
% %Keypoint.Sign(RemoveTooSmall,:) = [];
% Keypoint.Count = length(Keypoint.LocationIndex);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove Keypoints that are spatially close together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KeypointTree = KDTreeSearcher(Keypoint.Location, 'Distance', 'euclidean');

Ind = sub2ind(size(DoG), Keypoint.LocationIndex, Keypoint.Level);
KpDoG = DoG(Ind);

[KpDoGSort, KpDoGInd] = sort(abs(KpDoG),'descend');

RemoveInd = [];
RemoveIndAll = [];

KpIndRemain = ones(Keypoint.Count,1);
Discontinue = 0;

% WaitBar = waitbar(0, sprintf('Checking Keypoint %i of %i', 0, Keypoint.Count));


for i = 1 : Keypoint.Count

    CurrentPoint = Keypoint.Location(KpDoGInd(i),:);
    CurrentScale = Keypoint.Scale(KpDoGInd(i));
    CurrentDoG = KpDoGSort(i);
    CurrentInd = KpDoGInd(i);
    srange = t_range * CurrentScale;
%     srange = t_range * PointCloud.ResolutionLocal(Keypoint.Location(KpDoGInd(i)));
    
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
                        
                        % Marr-Hildreth 3D Edge Detection
                        if strcmpi(DoGNormalize, 'AbsDoG') && strcmpi(CompareMethod, '<')
                            if abs(NextDoG) > abs(CurrentDoG)
                                
                                RemoveInd = [RemoveInd, Neigh(j)];
                                
                            elseif abs(NextDoG) <= abs(CurrentDoG)
                                
                                RemoveInd = [RemoveInd, KpDoGInd(i)];
                                Discontinue = 1;
                                
                            else
                                warning('NMS needs another criterion.')
                            end
                            
                        % DoG local extrema
                        elseif strcmpi(DoGNormalize, 'NLoG') && strcmpi(CompareMethod, '<>')
                            
                            if abs(NextDoG) > abs(CurrentDoG)
                                
                                RemoveInd = [RemoveInd, KpDoGInd(i)];
                                Discontinue = 1;
                                
                            elseif abs(NextDoG) <= abs(CurrentDoG)
                                
                                RemoveInd = [RemoveInd, Neigh(j)];
                                
                                
                            else
                                warning('NMS needs another criterion.')
                            end
                            
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
    
%     waitbar(i/Keypoint.Count, WaitBar, sprintf('Checking Keypoint %i of %i', i, Keypoint.Count));
    
end

% waitbar(i/Keypoint.Count, WaitBar, sprintf('Checking Keypoints Complete'));
% close(WaitBar)


NMSKeypoint = Keypoint;


NMSKeypoint = rmfield(NMSKeypoint, 'Count');
FNames = fieldnames(NMSKeypoint);
for i = 1 : length(FNames)
    NMSKeypoint.(FNames{i})(RemoveIndAll,:) = [];
end
NMSKeypoint.Count = length(NMSKeypoint.LocationIndex);


% NMSKeypoint.Location(RemoveIndAll,:) = [];
% NMSKeypoint.Normal(RemoveIndAll,:) = [];
% NMSKeypoint.LocationIndex(RemoveIndAll,:) = [];
% NMSKeypoint.Scale(RemoveIndAll,:) = [];
% NMSKeypoint.Level(RemoveIndAll,:) = [];
% %NMSKeypoint.Sign(RemoveIndAll,:) = [];
% NMSKeypoint.Count = length(NMSKeypoint.LocationIndex);


end







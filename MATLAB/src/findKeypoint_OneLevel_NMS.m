% Andrew Rhodes
% WVU
% Jan. 2018
%
% findKeyPoints finds the keypoint locations and their scale on 2D images
% or 3D meshes.
%
% input[DoG]: Difference of Gaussian of signals over which to search
% input[ScaleParameter]: list of scale parameters accompanying the DoG
% input[Neighbors]: comparison neighbors for determining local extrema
% input[NumPoints]: the number of points that are being searched
% input[varargin]: dimension of data, either 2 or 3. Default is 3
%
% output[KeyPoint]: structure with elements .Location and .Scale

% (1) Compare to first ring neighbors
% (2) Require to be top NMSPercent
% (3) NMS within 3 sigma_i


function [Keypoint, NMSKeypoint] = findKeypoint_OneLevel_NMS(DoG, PointCloud, ScaleParameter, Neighbors, Method, CompareMethod, DoGNormalize, NMSPercent, t_range, t_scale)



if strcmpi(Method, 'Old')
    
    [NumVertices, NumDoGLevels] = size(DoG);

    if ~iscell(Neighbors)
        error('findKeypoint_OneLevel_NMS::Neighbors must be in cell format')
    end

    MaxNumberFeatures = round(NumDoGLevels * NumVertices / 100);
    
    
    % Initialize space for the keypoints
    
    
    Keypoint.Scale = zeros(MaxNumberFeatures,1);
    Keypoint.Location = zeros(MaxNumberFeatures,1);
    Keypoint.Level = zeros(MaxNumberFeatures,1);
    Keypoint.Response = zeros(MaxNumberFeatures,1);
    
    NumFeatures = 0;
    
    isFeaturePoint = 0;
    
    WaitBar = waitbar(0, sprintf('Checking Vertex %i of %i', 0, NumVertices));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% (1) Compare to first ring neighbors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : NumVertices
        CurrentNeighbors = Neighbors{i,1};
        
        for j = 2 : NumDoGLevels - 2
            
            CurrentValue = DoG(i,j);
            SurroundingValuesCurrent = DoG(CurrentNeighbors,j);
            
            if strcmpi(CompareMethod, '>')
                
                if all(CurrentValue > SurroundingValuesCurrent)
                    isFeaturePoint = 1;
                end
                
            elseif strcmpi(CompareMethod, '<>') || strcmpi(CompareMethod, '><')
                
                if all(CurrentValue > SurroundingValuesCurrent)
                    isFeaturePoint = 1;
                elseif all(CurrentValue < SurroundingValuesCurrent)                    
                    isFeaturePoint = 1;
                end
                
            elseif strcmpi(CompareMethod, '<')

                if all(CurrentValue < SurroundingValuesCurrent)
                    isFeaturePoint = 1;
                end
                
            else 
                error('findKeypoint_OneLevel_NMS:: incorrecte input for CompareMethod')
            end
            
            
            if isFeaturePoint
                NumFeatures = NumFeatures + 1;
                Keypoint.Scale(NumFeatures, 1) = ScaleParameter(j) + PointCloud.ResolutionLocal(i,1);
                Keypoint.Location(NumFeatures, 1) = i;
                Keypoint.Level(NumFeatures, 1) = j;
                Keypoint.Response(NumFeatures, 1) = CurrentValue;
                isFeaturePoint = 0;
            end
            
        end
        waitbar(i/NumVertices, WaitBar, sprintf('Checking Vertex %i of %i', i, NumVertices));
        
    end
    
    waitbar(i/NumVertices, WaitBar, sprintf('Checking Vertex for Keypoint Complete'));
    close(WaitBar)
    
    
    % Remove empty elements from initilization
    ZeroLogic = (Keypoint.Location == 0);
    
    Keypoint.Scale(ZeroLogic,:) = [];
    Keypoint.Location(ZeroLogic,:) = [];
    Keypoint.Level(ZeroLogic,:) = [];
    Keypoint.Response(ZeroLogic,:) = [];
    
    elseif strcmpi(Method, 'New')
        error('findKeypoint_OneLevel_NMS::Neighbors must be in cell format')
        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Require to be top NMSPercent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ResponseValue, ResponseOrder] = sort( Keypoint.Response, 'descend');

PercentPositions = round(NMSPercent * length(Keypoint.Scale)) + 1 : length(Keypoint.Scale);

% Remove the bottom (1-NMSPercent) responses and keypoints
ResponseValue(PercentPositions) = [];
ResponseOrder(PercentPositions) = [];


Keypoint.Scale = Keypoint.Scale(ResponseOrder);
Keypoint.Location = Keypoint.Location(ResponseOrder);
Keypoint.Level = Keypoint.Level(ResponseOrder);
Keypoint.Response = Keypoint.Response(ResponseOrder);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) NMS within 3 sigma_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



NumKeypoints = length(Keypoint.Location);

KeypointTree = KDTreeSearcher(PointCloud.Location(Keypoint.Location,:), 'Distance', 'euclidean');

Ind = sub2ind(size(DoG), Keypoint.Location, Keypoint.Level);
KpDoG = DoG(Ind);

[KpDoGSort, KpDoGInd] = sort(abs(KpDoG),'descend');

RemoveInd = [];
RemoveIndAll = [];

KpIndRemain = ones(NumKeypoints,1);
Discontinue = 0;

WaitBar = waitbar(0, sprintf('Checking Keypoint %i of %i', 0, NumKeypoints));


for i = 1 : NumKeypoints

    CurrentPoint = PointCloud.Location(Keypoint.Location(KpDoGInd(i)),:);
    CurrentScale = Keypoint.Scale(KpDoGInd(i));
    CurrentDoG = KpDoGSort(i);
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
                        
                        % Marr-Hildreth 3D Edge Detection
                        if strcmpi(DoGNormalize, 'AbsDoG') && strcmpi(CompareMethod, '<')
                            if abs(NextDoG) > abs(CurrentDoG)
                                
                                RemoveInd = [RemoveInd, Neigh(j)];
                                
                            elseif abs(NextDoG) <= abs(CurrentDoG)
                                
                                RemoveInd = [RemoveInd, KpDoGInd(i)];
                                Discontinue = 1;
                                
                            else
                                warning('findKeypoint_OneLevel_NMS::NMS needs another criterion.')
                            end
                            
                        % DoG local extrema
                        elseif strcmpi(DoGNormalize, 'DoG') && strcmpi(CompareMethod, '<>')
                            
                            if abs(NextDoG) > abs(CurrentDoG)
                                
                                RemoveInd = [RemoveInd, KpDoGInd(i)];
                                Discontinue = 1;
                                
                            elseif abs(NextDoG) <= abs(CurrentDoG)
                                
                                RemoveInd = [RemoveInd, Neigh(j)];
                                
                                
                            else
                                warning('findKeypoint_OneLevel_NMS::NMS needs another criterion.')
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
    
    waitbar(i/NumKeypoints, WaitBar, sprintf('Checking Keypoint %i of %i', i, NumKeypoints));
    
        
        
end

waitbar(i/NumKeypoints, WaitBar, sprintf('Checking Keypoints Complete'));
close(WaitBar)



NMSKeypoint = Keypoint;

NMSKeypoint.Location(RemoveIndAll) = [];
NMSKeypoint.Scale(RemoveIndAll) = [];
NMSKeypoint.Level(RemoveIndAll) = [];
NMSKeypoint.Response(RemoveIndAll) = [];

end






























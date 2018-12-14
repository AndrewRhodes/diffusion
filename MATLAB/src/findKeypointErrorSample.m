





function [Error, NoMatch] = findKeypointErrorSample...
    (PointCloudOriginal, KeypointOriginal, NumStepsOriginal, PointCloud,...
                    KeypointIter, NumSteps, level_min, t_scale, t_range, tauFraction)
                
    alpha = 1;

if level_min
    ScaleLogic = KeypointOriginal.Scale < ( level_min * PointCloudOriginal.ResolutionLocal(KeypointOriginal.LocationIndex) );
    ScaleLogic2 = KeypointOriginal.Scale < PointCloud.Resolution;
%     ScaleLogic2 = KeypointOriginal.Scale < PointCloud.Resolution + level_min * PointCloud.Resolution;
    
    ScaleLogic = ScaleLogic | ScaleLogic2;
    if isfield(KeypointOriginal, 'Count')
        KeypointOriginal = rmfield(KeypointOriginal, 'Count');
    end
    FNames = fieldnames(KeypointOriginal);
    for jj = 1 : length(FNames)
        KeypointOriginal.(FNames{jj})(ScaleLogic,:) = [];
    end
    KeypointOriginal.Count = length(KeypointOriginal.LocationIndex);
    KeypointOriginal.Count
    
    
    ScaleLogic = KeypointIter.Scale < ( level_min * PointCloud.ResolutionLocal(KeypointIter.LocationIndex) );
    if isfield(KeypointOriginal, 'Count')
        KeypointIter = rmfield(KeypointIter, 'Count');
    end
    FNames = fieldnames(KeypointIter);
    for jj = 1 : length(FNames)
        KeypointIter.(FNames{jj})(ScaleLogic,:) = [];
    end
    KeypointIter.Count = length(KeypointIter.LocationIndex);
    KeypointIter.Count
    
    
end    






Error.Scale = 0;
Error.Distance = 0;
Error.Count = 0;
Error.RelativeRepeat = 0;
Error.ScaleRepeat = 0;

NoMatch = 0;



KeypointIterTree = KDTreeSearcher(KeypointIter.Location, 'Distance', 'euclidean');

c = 0;
Matched = zeros(KeypointIter.Count,1);

for j = 1 : KeypointOriginal.Count

    CurrentPoint = KeypointOriginal.Location(j,:);
    CurrentScale = KeypointOriginal.Scale(j,1);

	srange = t_range * CurrentScale;
% 	srange = 2 * PointCloud.Resolution;
        
    [Neigh, Dist] = rangesearch(KeypointIterTree, CurrentPoint, srange);
        
    Neigh = Neigh{:};
    Dist = Dist{:};
               
    % Remove from consideration the original keypoints that have
    % already been matched.
    Remove = [];
    for k = 1 : length(Neigh)
        if Matched(Neigh)
            Remove = [Remove; k];
        end
    end
    if ~isempty(Neigh)
        Neigh(Remove) = [];
        Dist(Remove) = [];
    end    
        
    if isempty(Neigh)
        NoMatch = NoMatch + 1;  
        
    else
        ScaleRatio = min(KeypointIter.Scale(Neigh), CurrentScale) ./ max(KeypointIter.Scale(Neigh), CurrentScale) ;
        ScaleLogic = ScaleRatio >= t_scale;

        if nnz(ScaleLogic) == 0
            NoMatch = NoMatch + 1;
        else
            [MaxVal, MaxLoc] = max(ScaleRatio);
            c = c + 1;
            Error.Scale(c,1) = findSphereRepeat(KeypointIter.Location(Neigh(MaxLoc),:), KeypointIter.Scale(Neigh(MaxLoc)), CurrentPoint, CurrentScale);
            Error.Distance(c,1) = Dist(MaxLoc);
            Error.Count = Error.Count + 1;
            Matched(Neigh(MaxLoc)) = 1;        
        end
    end
        
    if nnz(Matched) == KeypointIter.Count
        warning('All Keypoint Matched Early')
        break;
    end
        
end


Error.RelativeRepeat = Error.Count / KeypointOriginal.Count;
Error.ScaleRepeat = sum(Error.Scale) / Error.Count;
    
    


end




function Vratio = findSphereRepeat(p1, r1, p2, r2)
% Vratio = Vintersect / Vunion

rmin = min(r1,r2);
rmax = max(r1,r2);

d = sqrt(sum((p1-p2).^2));

if d > r1 + r2
    Vratio = 0;
    return;
elseif r1 >= d + r2
    %     Sphere r2 is completely inside Sphere r1
    Vintersect = findSphereVolume(r2);
    Vunion = findSphereVolume(r1);
elseif r2 >= d + r1
    %     Sphere r1 is completely inside Sphere r2
    Vintersect = findSphereVolume(r1);
    Vunion = findSphereVolume(r2);
else % Spheres intersect
   
    x = (d^2 - rmin^2 + rmax^2) / (2*d);
    h1 = rmax - x;
    h2 = rmin - d + x;
    
    V1cap = findCapVolume(rmax, h1);
    V2cap = findCapVolume(rmin, h2);
    
    Vintersect = V1cap + V2cap;
    
    V1 = findSphereVolume(r1);
    V2 = findSphereVolume(r2);
    
    Vunion = V1 - V1cap + V2 - V2cap;
end

Vratio = Vintersect / Vunion;

end




function V = findCapVolume(radius, height)

V = (1/3) * pi * height^2 * (3*radius-height);


end




function V = findSphereVolume(radius)

V = (4/3)*pi*radius^3;

end







function [Error, NoMatch] = findKeypointErrorNew(PointCloud, FileLocation, NumIter, UseNMS, t_scale, level_min, t_range)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model File Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/bunny/');
if UseNMS
    FileName = strcat('NMSKeypoint','.mat');
    load(fullfile(FileLocation, FileName),'NMSKeypoint')
    KeypointOriginal = NMSKeypoint;
else
    FileName = strcat('Keypoint','.mat');
    load(fullfile(FileLocation, FileName),'Keypoint')
    KeypointOriginal = Keypoint;
end

if level_min
    ScaleLogic = KeypointOriginal.Level < level_min;
    KeypointOriginal.Scale(ScaleLogic) = [];
    KeypointOriginal.Location(ScaleLogic) = [];
    KeypointOriginal.Level(ScaleLogic) = [];
end

NumKeypointOriginal = length(KeypointOriginal.Scale)


KeypointTree = KDTreeSearcher(PointCloud.Location(KeypointOriginal.Location,:), 'Distance', 'euclidean');



Error.Scale = cell(NumIter,1);
Error.Distance = cell(NumIter,1);
Error.Count = zeros(NumIter,1);
Error.RelativeRepeat = zeros(NumIter,1);
Error.ScaleRepeat = zeros(NumIter,1);
 
 

NoMatch = zeros(NumIter,1);
 

for i = 1 : NumIter
    
    if UseNMS
        FileName = strcat('NMSKeypoint','_Iter',num2str(i),'.mat');
        load(fullfile(FileLocation, FileName),'NMSKeypoint')
        KeypointIter = NMSKeypoint;
    else
        FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
        load(fullfile(FileLocation, FileName),'Keypoint')
        KeypointIter = Keypoint;
    end
    
    
    if level_min
        ScaleLogic = KeypointIter.Level < level_min;
        KeypointIter.Scale(ScaleLogic) = [];
        KeypointIter.Location(ScaleLogic) = [];
        KeypointIter.Level(ScaleLogic) = [];
    end
        
    
    NumKeypointsIter = length(KeypointIter.Location);
    c = 0;
    Matched = zeros(NumKeypointOriginal,1);
    
    for j = 1 : NumKeypointsIter
        
        CurrentPoint = PointCloud.Location(KeypointIter.Location(j),:);
        CurrentScale = KeypointIter.Scale(j);
        srange = t_range * PointCloud.ResolutionLocal(KeypointIter.Location(j));
        
        % [Neigh, Dist] = rangesearch(KeypointTree, CurrentPoint, PointCloud.Resolution);
        [Neigh, Dist] = rangesearch(KeypointTree, CurrentPoint, srange);
        % [Idx, Dist] = knnsearch(KeypointTree, CurrentPoint, 'K', 1, 'IncludeTies', true, 'Distance', 'euclidean');

        
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
            
            NoMatch(i,1) = NoMatch(i,1) + 1;
%             continue;
            
        elseif length(Neigh) == 1
            
            ScaleRatio = min(KeypointOriginal.Scale(Neigh), CurrentScale) / max(KeypointOriginal.Scale(Neigh), CurrentScale) ;
            
            if ScaleRatio >= t_scale
                c = c + 1;
                Error.Scale{i}(c,1) = findSphereRepeat(PointCloud.Location(KeypointOriginal.Location(Neigh),:), KeypointOriginal.Scale(Neigh), CurrentPoint, CurrentScale);
                Error.Distance{i}(c,1) = Dist;
                Error.Count(i,1) = Error.Count(i,1) + 1;
%                 Match(i,1) = Match(i,1) + 1;
                Matched(Neigh) = 1;
            else
                NoMatch(i,1) = NoMatch(i,1) + 1;
            end
           
        elseif length(Neigh) > 1
            
            ScaleRatio = min(KeypointOriginal.Scale(Neigh), CurrentScale) ./ max(KeypointOriginal.Scale(Neigh), CurrentScale) ;
            ScaleLogic = ScaleRatio >= t_scale;
            
            if nnz(ScaleLogic) == 0
                NoMatch(i,1) = NoMatch(i,1) + 1;
            else
                [MaxVal, MaxLoc] = max(ScaleRatio);
                c = c + 1;
                Error.Scale{i}(c,1) = findSphereRepeat(PointCloud.Location(KeypointOriginal.Location(Neigh(MaxLoc)),:), KeypointOriginal.Scale(Neigh(MaxLoc)), CurrentPoint, CurrentScale);
                Error.Distance{i}(c,1) = Dist(MaxLoc);
                Error.Count(i,1) = Error.Count(i,1) + 1;
%                 Match(i,1) = Match(i,1) + 1;
                Matched(Neigh(MaxLoc)) = 1;
            end
            
        end
        
    if nnz(Matched) == NumKeypointOriginal
        error('Ended Early')
    end
            
  
    end
    
    Error.RelativeRepeat(i,1) = Error.Count(i,1) / NumKeypointOriginal;
    Error.ScaleRepeat(i,1) = sum(Error.Scale{i,1}) / Error.Count(i,1);
    
   
    
end





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

































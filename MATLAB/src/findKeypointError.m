


function [Error, Match, NoMatch, MultipleMatch] = findKeypointError(PointCloud, FileLocation, NumIter, UseNMS, Use4D, t_scale, level_min)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UseNMS = 0;
% t_scale = 0.9;


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

NumKeypointOriginal = length(KeypointOriginal.Scale);


if Use4D
    KeypointTree = KDTreeSearcher([PointCloud.Location(KeypointOriginal.Location,:),KeypointOriginal.Scale], 'Distance', 'euclidean');    
else
    KeypointTree = KDTreeSearcher(PointCloud.Location(KeypointOriginal.Location,:), 'Distance', 'euclidean');
end


Error.Scale = cell(NumIter,1);
Error.Distance = cell(NumIter,1);
Error.Count = zeros(NumIter,1);
Error.RelativeRepeat = zeros(NumIter,1);
Error.ScaleRepeat = zeros(NumIter,1);
 
 

NoMatch = zeros(NumIter,1);
Match = zeros(NumIter,1);
MultipleMatch = zeros(NumIter,1);
 

for i = 1 : NumIter
    
    if UseNMS
        FileName = strcat('NMSKeypoint','_Iter',num2str(i),'.mat');
        load(fullfile(FileLocation, FileName),'NMSKeypoint')
        Keypoint = NMSKeypoint;
    else
        FileName = strcat('Keypoint','_Iter',num2str(i),'.mat');
        load(fullfile(FileLocation, FileName),'Keypoint')
    end
    
    if level_min
        ScaleLogic = Keypoint.Level < level_min;
        Keypoint.Scale(ScaleLogic) = [];
        Keypoint.Location(ScaleLogic) = [];
        Keypoint.Level(ScaleLogic) = [];
    end
        
    
    NumKeypoints = length(Keypoint.Location);
    c = 0;
    for j = 1 : NumKeypoints
        
        CurrentPoint = PointCloud.Location(Keypoint.Location(j),:);
        CurrentScale = Keypoint.Scale(j);
        
        if Use4D
%             [Neigh, Dist] = rangesearch(KeypointTree, [CurrentPoint,CurrentScale], PointCloud.Resolution);
            [Neigh, Dist] = rangesearch(KeypointTree, [CurrentPoint,CurrentScale], CurrentScale);
%             [Idx, Dist] = knnsearch(KeypointTree, [CurrentPoint,CurrentScale], 'K', 1, 'IncludeTies', true, 'Distance', 'euclidean');
        else
%             [Neigh, Dist] = rangesearch(KeypointTree, CurrentPoint, PointCloud.Resolution);
            [Neigh, Dist] = rangesearch(KeypointTree, CurrentPoint, CurrentScale);
%             [Idx, Dist] = knnsearch(KeypointTree, CurrentPoint, 'K', 1, 'IncludeTies', true, 'Distance', 'euclidean');
        end
        
        Neigh = Neigh{:};
        Dist = Dist{:};
        MoreThanOne = 0;
        OnlyOne = 0;
        
        if isempty(Dist)
%             warning('No Match on Iter %d, Keypoint %d.', i, j )
            NoMatch(i,1) = NoMatch(i,1) + 1;
            continue;
        elseif length(Dist) == 1
            OnlyOne = 1;
            ScaleRatio = (min(KeypointOriginal.Scale(Neigh), CurrentScale) ./ max(KeypointOriginal.Scale(Neigh), CurrentScale));
            ScaleLogic = ScaleRatio >= t_scale;            
        elseif length(Dist) >= 1
            MoreThanOne = 1;
            ScaleRatio = (min(KeypointOriginal.Scale(Neigh), CurrentScale) ./ max(KeypointOriginal.Scale(Neigh), CurrentScale));
            ScaleLogic = ScaleRatio >= t_scale;            
        else
            warning('Didn''t think of something')
        end
        
        
        
        
        if MoreThanOne && nnz(ScaleLogic) > 1
%             warning('Too many matches on Iter %d, Keypoint %d.', i, j)
%             c = c + 1;
            MultipleMatch(i,1) = MultipleMatch(i,1) + 1;
%             [maxval, maxind] = max(ScaleRatio);
%             Error.Scale{i}(c,1) = findSphereRepeat(PointCloud.Location(KeypointOriginal.Location(Neigh(maxind)),:), KeypointOriginal.Scale(Neigh(maxind)), CurrentPoint, CurrentScale);
%             Error.Distance{i}(c,1) = Dist(maxind);
%             Error.Count(i,1) = Error.Count(i,1) + 1;
        elseif MoreThanOne && nnz(ScaleLogic) == 0
%             warning('No Match on Iter %d, Keypoint %d.', i, j )
            NoMatch(i,1) = NoMatch(i,1) + 1;
        elseif MoreThanOne && nnz(ScaleLogic) == 1
            c = c + 1;
            Error.Scale{i}(c,1) = findSphereRepeat(PointCloud.Location(KeypointOriginal.Location(Neigh(ScaleLogic)),:), KeypointOriginal.Scale(Neigh(ScaleLogic)), CurrentPoint, CurrentScale);
            Error.Distance{i}(c,1) = Dist(ScaleLogic);
            Error.Count(i,1) = Error.Count(i,1) + 1;
            Match(i,1) = Match(i,1) + 1;
        elseif OnlyOne && nnz(ScaleLogic) == 1
            c = c + 1;
%             findSphereIntersect(PointCloud.Location(KeypointOriginal.Location(Neigh),:), KeypointOriginal.Scale(Neigh), CurrentPoint, CurrentScale)
            Error.Scale{i}(c,1) = findSphereRepeat(PointCloud.Location(KeypointOriginal.Location(Neigh),:), KeypointOriginal.Scale(Neigh), CurrentPoint, CurrentScale);
            Error.Distance{i}(c,1) = Dist;
            Error.Count(i,1) = Error.Count(i,1) + 1;
            Match(i,1) = Match(i,1) + 1;
        else % OnlyOne && nnz(ScaleLogic) > 1
%             warning('No Match on Iter %d, Keypoint %d.', i, j )
%             MoreThanOne
%             OnlyOne
%             ScaleLogic
%             (min(KeypointOriginal.Scale(Neigh), CurrentScale) ./ max(KeypointOriginal.Scale(Neigh), CurrentScale))
            NoMatch(i,1) = NoMatch(i,1) + 1;
            continue;
        end            
            
  
    end
    
    Error.RelativeRepeat(i,1) = Error.Count(i,1) / NumKeypointOriginal;
    Error.ScaleRepeat(i,1) = sum(Error.Scale{i,1}) / Error.Count(i,1);
%     Error.ScaleRepeat(i,1) = sum(Error.Scale{i}) / NumKeypointOriginal;
    
   
    
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

































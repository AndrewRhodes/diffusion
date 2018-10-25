




function [Error, NoMatch] = findKeypointErrorNewer(PointCloud, FileLocation, NoiseFileLocation, NumIter, t_scale, level_min, t_range)


load(strcat(FileLocation, 'Keypoint.mat'),'Keypoint')
load(strcat(FileLocation, 'NewKeypoint.mat'),'NewKeypoint')
load(strcat(FileLocation, 'Zeropoint.mat'),'Zeropoint')

NewKeypointTrue = NewKeypoint;

if level_min
    ScaleLogic = (NewKeypointTrue.Scale ./ PointCloud.ResolutionLocal(NewKeypointTrue.Location)) < level_min;
% 	ScaleLogic = ((NewKeypointTrue.Scale - PointCloud.ResolutionLocal(NewKeypointTrue.Location)) ./ PointCloud.Resolution) < level_min;

    NewKeypointTrue.Scale(ScaleLogic) = [];
    NewKeypointTrue.Location(ScaleLogic) = [];
    NewKeypointTrue.Level(ScaleLogic) = [];
    NewKeypointTrue.Count = length(NewKeypointTrue.Level);
end

[a,b] = unique(NewKeypointTrue.Location);
NewKeypointTrue.Scale = NewKeypointTrue.Scale(b);
NewKeypointTrue.Location = NewKeypointTrue.Location(b);
NewKeypointTrue.Level = NewKeypointTrue.Level(b);
NewKeypointTrue.Count = length(NewKeypointTrue.Level);


Tree = KDTreeSearcher(PointCloud.Location(NewKeypointTrue.Location,:), 'Distance', 'Euclidean');


Error.Scale = cell(NumIter,1);
Error.Distance = cell(NumIter,1);
Error.Count = zeros(NumIter,1);
Error.RelativeRepeat = zeros(NumIter,1);
Error.ScaleRepeat = zeros(NumIter,1);


NoMatch = zeros(NumIter,1);


for i = 1 : NumIter
    
    load(strcat(NoiseFileLocation, 'Keypoint_Iter',num2str(i),'.mat'),'Keypoint')
    load(strcat(NoiseFileLocation, 'NewKeypoint_Iter',num2str(i),'.mat'),'NewKeypoint')
    load(strcat(NoiseFileLocation, 'Zeropoint_Iter',num2str(i),'.mat'),'Zeropoint')
    
    if level_min
        ScaleLogic = (NewKeypoint.Scale ./ PointCloud.ResolutionLocal(NewKeypoint.Location)) < level_min;
%         ScaleLogic = ((NewKeypoint.Scale - PointCloud.ResolutionLocal(NewKeypoint.Location)) ./ PointCloud.Resolution) < level_min;
        
        NewKeypoint.Scale(ScaleLogic) = [];
        NewKeypoint.Location(ScaleLogic) = [];
        NewKeypoint.Level(ScaleLogic) = [];
        NewKeypoint.Count = length(NewKeypoint.Level);
    end
    [a,b] = unique(NewKeypoint.Location);
    NewKeypoint.Scale = NewKeypoint.Scale(b);
    NewKeypoint.Location = NewKeypoint.Location(b);
    NewKeypoint.Level = NewKeypoint.Level(b);
    NewKeypoint.Count = length(NewKeypoint.Level);
    
    
    c = 0;
    Matched = zeros(NewKeypointTrue.Count,1);
    
   for j = 1 : NewKeypoint.Count
       
       CurrentPoint = PointCloud.Location(NewKeypoint.Location(j),:);
       CurrentScale = NewKeypoint.Scale(j); % Absolute scale_i
       
       srange = t_range * CurrentScale;       
%        srange = t_range * PointCloud.ResolutionLocal(NewKeypoint.Location(j));
       
       [Neigh, Dist] = rangesearch(Tree, CurrentPoint, srange);
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
            
        elseif length(Neigh) == 1
            
            ScaleRatio = min(NewKeypointTrue.Scale(Neigh), CurrentScale) / max(NewKeypointTrue.Scale(Neigh), CurrentScale) ;
            
            if ScaleRatio >= t_scale
                c = c + 1;
                Error.Scale{i}(c,1) = findSphereRepeat(PointCloud.Location(NewKeypointTrue.Location(Neigh),:), NewKeypointTrue.Scale(Neigh), CurrentPoint, CurrentScale);
                Error.Distance{i}(c,1) = Dist;
                Error.Count(i,1) = Error.Count(i,1) + 1;
                Matched(Neigh) = 1;
            else
                NoMatch(i,1) = NoMatch(i,1) + 1;
            end
            
        elseif length(Neigh) > 1
            
            ScaleRatio = min(NewKeypointTrue.Scale(Neigh), CurrentScale) ./ max(NewKeypointTrue.Scale(Neigh), CurrentScale) ;
            ScaleLogic = ScaleRatio >= t_scale;
            
            if nnz(ScaleLogic) == 0
                NoMatch(i,1) = NoMatch(i,1) + 1;
            else
                [MaxVal, MaxLoc] = max(ScaleRatio);
                c = c + 1;
                Error.Scale{i}(c,1) = findSphereRepeat(PointCloud.Location(NewKeypointTrue.Location(Neigh(MaxLoc)),:), NewKeypointTrue.Scale(Neigh(MaxLoc)), CurrentPoint, CurrentScale);
                Error.Distance{i}(c,1) = Dist(MaxLoc);
                Error.Count(i,1) = Error.Count(i,1) + 1;
                Matched(Neigh(MaxLoc)) = 1;
            end
            
        end
         
        
        if nnz(Matched) == NewKeypointTrue.Count
            error('Ended Early')
        end
        
        
   end
   
	Error.RelativeRepeat(i,1) = Error.Count(i,1) / NewKeypointTrue.Count;
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

























function [Error, NoMatch] = keypointErrorTest(PointCloud, FileLocation, NoiseFileLocation, NumIter, UseNMS, t_scale, level_min, t_range, NoisyPCLocation, NoisyPCNeighbors)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model File Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/bunny/');

if UseNMS
    load(strcat(FileLocation, 'NMSKeypoint.mat'),'NMSKeypoint')
    KeypointOriginal = NMSKeypoint;
else
    load(strcat(FileLocation, 'Keypoint.mat'),'Keypoint')
    KeypointOriginal = Keypoint;
end

ScaleParameter = findScaleParameter(PointCloud.Resolution/12, 1/2, 3000, 'Laplacian', 'Natural');
KeypointOriginal.Scale = ScaleParameter(KeypointOriginal.Level) + PointCloud.ResolutionLocal(KeypointOriginal.LocationIndex);


if level_min
    ScaleLogic = KeypointOriginal.Scale < ( level_min * PointCloud.ResolutionLocal(KeypointOriginal.LocationIndex) );
    
    KeypointOriginal = rmfield(KeypointOriginal, 'Count');
    FNames = fieldnames(KeypointOriginal);
    for jj = 1 : length(FNames)
        KeypointOriginal.(FNames{jj})(ScaleLogic,:) = [];
    end
    KeypointOriginal.Count = length(KeypointOriginal.LocationIndex);
%     
%     ScaleLogic = KeypointOriginal.Level > 2000;
%     KeypointOriginal = rmfield(KeypointOriginal, 'Count');
%     FNames = fieldnames(KeypointOriginal);
%     for jj = 1 : length(FNames)
%         KeypointOriginal.(FNames{jj})(ScaleLogic,:) = [];
%     end
%     KeypointOriginal.Count = length(KeypointOriginal.LocationIndex);

    % % % % %     ScaleLogic = ((KeypointOriginal.Scale - PointCloud.ResolutionLocal(KeypointOriginal.Location)) ./ PointCloud.Resolution) < level_min;
    % % % % %     ScaleLogic = (KeypointOriginal.Scale ./ PointCloud.Resolution) < 4;
    %     ScaleLogic = KeypointOriginal.Level > level_min;
%     KeypointOriginal.Scale(ScaleLogic) = [];
%     KeypointOriginal.Location(ScaleLogic) = [];
%     KeypointOriginal.Level(ScaleLogic) = [];
end

% KeypointOriginal.Count

% KeypointTree = KDTreeSearcher(KeypointOriginal.Location, 'Distance', 'euclidean');



Error.Scale = cell(NumIter,1);
Error.Distance = cell(NumIter,1);
Error.Count = zeros(NumIter,1);
Error.RelativeRepeat = zeros(NumIter,1);
Error.ScaleRepeat = zeros(NumIter,1);



NoMatch = zeros(NumIter,1);


for i = 1 : NumIter
    
    [PointCloudNoisy.Location, PointCloudNoisy.Face, PointCloudNoisy.Normal, PointCloudNoisy.Signal]...
                = read_ply_all_elements( strcat(NoisyPCLocation,'_iter',num2str(i),'.ply') );
    PointCloudNoisy.LocationCount = size(PointCloudNoisy.Location,1);
    PointCloudNoisy.FaceCount = size(PointCloudNoisy.Face, 1);
    
    load( strcat(NoisyPCNeighbors,'_iter',num2str(i),'.mat'), 'Neighbors')
    PointCloudNoisy = findLocalResolution(PointCloudNoisy, Neighbors.Connect);
    PointCloudNoisy = findMeshResolution(PointCloudNoisy, 'Model');
    ScaleParameter = findScaleParameter(PointCloudNoisy.Resolution/12, 1/2, 3000, 'Laplacian', 'Natural');
    
    
    if UseNMS
        load(strcat(NoiseFileLocation, 'NMSKeypoint','_Iter',num2str(i),'.mat'),'NMSKeypoint')
        KeypointIter = NMSKeypoint;
    else
        load(strcat(NoiseFileLocation, 'Keypoint','_Iter',num2str(i),'.mat'),'Keypoint')
        KeypointIter = Keypoint;
    end
    
    KeypointIter.Scale = ScaleParameter(KeypointIter.Level) + PointCloudNoisy.ResolutionLocal(KeypointIter.LocationIndex);
    
%     KeypointIter.Count
    
    if level_min
        ScaleLogic = KeypointIter.Scale < ( level_min * PointCloudNoisy.ResolutionLocal(KeypointIter.LocationIndex));
        
        KeypointIter = rmfield(KeypointIter, 'Count');
        FNames = fieldnames(KeypointIter);
        for jj = 1 : length(FNames)
            KeypointIter.(FNames{jj})(ScaleLogic,:) = [];
        end
        KeypointIter.Count = length(KeypointIter.LocationIndex);

        
%         ScaleLogic = KeypointIter.Level > 2000;
%         KeypointIter = rmfield(KeypointIter, 'Count');
%         FNames = fieldnames(KeypointIter);
%         for jj = 1 : length(FNames)
%             KeypointIter.(FNames{jj})(ScaleLogic,:) = [];
%         end
%         KeypointIter.Count = length(KeypointIter.LocationIndex);
    
    
        % % % % %         ScaleLogic = ((KeypointIter.Scale - PointCloud.ResolutionLocal(KeypointIter.Location)) ./ PointCloud.Resolution) < level_min;
        % % % % %         ScaleLogic = (KeypointIter.Scale ./ PointCloud.Resolution) < 4;
        %         ScaleLogic = KeypointIter.Level > level_min;
%         KeypointIter.Scale(ScaleLogic) = [];
%         KeypointIter.Location(ScaleLogic,:) = [];
%         KeypointIter.Level(ScaleLogic) = [];
    end
    
    KeypointIterTree = KDTreeSearcher(KeypointIter.Location, 'Distance', 'euclidean');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Replace noisy PC vertice with vertex locations of original PC
    KeypointIter.Location = PointCloud.Location(KeypointIter.LocationIndex,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    c = 0;
    Matched = zeros(KeypointIter.Count,1);
        
    for j = 1 : KeypointOriginal.Count
        CurrentPoint = KeypointOriginal.Location(j,:);
        CurrentScale = KeypointOriginal.Scale(j,1);
        srange = t_range * CurrentScale; % PointCloud.Resolution * 3; %
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
            NoMatch(i,1) = NoMatch(i,1) + 1;
        
        else
            ScaleRatio = min(KeypointIter.Scale(Neigh), CurrentScale) ./ max(KeypointIter.Scale(Neigh), CurrentScale) ;
            ScaleLogic = ScaleRatio >= t_scale;
            
            if nnz(ScaleLogic) == 0
                NoMatch(i,1) = NoMatch(i,1) + 1;
            else
                [MaxVal, MaxLoc] = max(ScaleRatio);
                c = c + 1;
                Error.Scale{i}(c,1) = findSphereRepeat(KeypointIter.Location(Neigh(MaxLoc),:), KeypointIter.Scale(Neigh(MaxLoc)), CurrentPoint, CurrentScale);
                Error.Distance{i}(c,1) = Dist(MaxLoc);
                Error.Count(i,1) = Error.Count(i,1) + 1;
                %                 Match(i,1) = Match(i,1) + 1;
                Matched(Neigh(MaxLoc)) = 1;
            end
            
        end
        
        
        if nnz(Matched) == KeypointIter.Count
            warning('All Keypoint Matched Early')
            break;
        end
        
    end
    
    Error.RelativeRepeat(i,1) = Error.Count(i,1) / KeypointOriginal.Count;
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

































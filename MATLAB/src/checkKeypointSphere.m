




function NewKeypoint = checkKeypointSphere(PointCloud, ScaleParameter, Keypoint, Zeropoint, t, Method)

NumLevels = length(ScaleParameter);

NumKeypoints = length(Keypoint.Location);

KeypointTree = KDTreeSearcher(PointCloud.Location, 'Distance', 'euclidean');

NumNewKeypoint = 0;
NewKeypoint.Level = zeros(NumKeypoints,1);
NewKeypoint.Location = zeros(NumKeypoints,1);
NewKeypoint.Scale = zeros(NumKeypoints,1);
NewKeypoint.Sign = zeros(NumKeypoints,1);
NewKeypoint.DoG = zeros(NumKeypoints,1);
NewKeypoint.Count = 0;

sp = 0;
ep = 0;

WaitBar = waitbar(0, sprintf('Checking Level %i of %i', 0, NumLevels));

if strcmpi(Method, 'Local')
       
    searcher = @(Tree, Point, Radius) rangesearch(Tree, Point, Radius);
    
    for i = 1 : NumLevels - 1
        
        KeypointsAtLevel = find(Keypoint.Level == i);
        NumKeypointsAtLevel = length(KeypointsAtLevel);
        
        ZeropointAtLevel = Zeropoint.Location(Zeropoint.Level == i);
        
        srange1 = ScaleParameter(i) + PointCloud.ResolutionLocal(Keypoint.Location(KeypointsAtLevel));
        srange2 = ScaleParameter(i+1) + PointCloud.ResolutionLocal(Keypoint.Location(KeypointsAtLevel));
        
        CurrentPoint = PointCloud.Location(Keypoint.Location(KeypointsAtLevel),:);
        
        Neigh1 = cellfun(searcher, repmat({KeypointTree}, NumKeypointsAtLevel,1), mat2cell(CurrentPoint,ones(NumKeypointsAtLevel,1), 3), num2cell(srange1), 'Uni', 0);
        Neigh2 = cellfun(searcher, repmat({KeypointTree}, NumKeypointsAtLevel,1), mat2cell(CurrentPoint,ones(NumKeypointsAtLevel,1), 3), num2cell(srange2), 'Uni', 0);
        
        if ~isempty(Neigh1)
            Neigh1 = [Neigh1{:}]';
        end
        if ~isempty(Neigh2)
            Neigh2 = [Neigh2{:}]';
        end

        ZeroCrossingCell = repmat( {ZeropointAtLevel}, NumKeypointsAtLevel, 1);
        
        ZeroSphere1 = cellfun(@ismember, Neigh1, ZeroCrossingCell, 'Uni', 0);
        ZeroSphere2 = cellfun(@ismember, Neigh2, ZeroCrossingCell, 'Uni', 0);
        
        nnzZeroSphere1 = cell2mat(cellfun(@nnz, ZeroSphere1, 'Uni', 0));
        nnzZeroSphere2 = cell2mat(cellfun(@nnz, ZeroSphere2, 'Uni', 0));
        
        Locations = find((nnzZeroSphere1 == 0) & (nnzZeroSphere2 >0));
        
        if ~isempty(Locations)
            sp = ep + 1;
            ep = sp + length(Locations) - 1;
            NewKeypoint.Level(sp:ep,1) = i;
            NewKeypoint.Location(sp:ep, 1) = Keypoint.Location(KeypointsAtLevel(Locations));
            NewKeypoint.Scale(sp:ep, 1) = srange1(Locations) ./ t;
            NewKeypoint.Count = NewKeypoint.Count + length(Locations);
            NewKeypoint.Sign(sp:ep, 1) = Keypoint.Sign(KeypointsAtLevel(Locations));
            NewKeypoint.DoG(sp:ep, 1) = Keypoint.DoG(KeypointsAtLevel(Locations));
        end
        
        waitbar(i/NumLevels, WaitBar, sprintf('Checking Level %i of %i', i, NumLevels));
        
    end
    
    
    
elseif strcmpi(Method, 'Global')
    
    
    for i = 1 : NumLevels - 1
        
        KeypointsAtLevel = find(Keypoint.Level == i);
        NumKeypointsAtLevel = length(KeypointsAtLevel);
        
        
        ZeropointAtLevel = Zeropoint.Location(Zeropoint.Level == i);
                
        srange1 = t * (ScaleParameter(i) + PointCloud.Resolution);
        srange2 = t * (ScaleParameter(i+1) + PointCloud.Resolution);
        
        CurrentPoint = PointCloud.Location(Keypoint.Location(KeypointsAtLevel),:);
        
        [Neigh1, ~] = rangesearch(KeypointTree, CurrentPoint, srange1);
        [Neigh2, ~] = rangesearch(KeypointTree, CurrentPoint, srange2);
        
        ZeroCrossingCell = repmat( {ZeropointAtLevel}, length(CurrentPoint), 1);
        
        ZeroSphere1 = cellfun(@ismember, Neigh1, ZeroCrossingCell, 'Uni', 0);
        ZeroSphere2 = cellfun(@ismember, Neigh2, ZeroCrossingCell, 'Uni', 0);
        
        nnzZeroSphere1 = cell2mat(cellfun(@nnz, ZeroSphere1, 'Uni', 0));
        nnzZeroSphere2 = cell2mat(cellfun(@nnz, ZeroSphere2, 'Uni', 0));
        
        Locations = find((nnzZeroSphere1 == 0) & (nnzZeroSphere2 >0));
        
        if ~isempty(Locations)
            sp = ep + 1;
            ep = sp + length(Locations) - 1;
            NewKeypoint.Level(sp:ep,1) = i;
            NewKeypoint.Location(sp:ep, 1) = Keypoint.Location(KeypointsAtLevel(Locations));
            NewKeypoint.Scale(sp:ep, 1) = srange1 / t;
            NewKeypoint.Count = NewKeypoint.Count + length(Locations);
            NewKeypoint.Sign(sp:ep, 1) = Keypoint.Sign(KeypointsAtLevel(Locations));
            NewKeypoint.DoG(sp:ep, 1) = Keypoint.DoG(KeypointsAtLevel(Locations));
        end
        
        waitbar(i/NumLevels, WaitBar, sprintf('Checking Level %i of %i', i, NumLevels));
        
    end
    
end

waitbar(i/NumLevels, WaitBar, sprintf('Checking Levels Complete'));
close(WaitBar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 1 : NumLevels - 1
%
%     KeypointsAtLevel = find(Keypoint.Level == i);
%
%     ZeropointAtLevel = Zeropoint.Location(Zeropoint.Level == i);
%
%     for j = 1 : length(KeypointsAtLevel)
%
%         CurrentPoint = PointCloud.Location(Keypoint.Location(KeypointsAtLevel(j)),:);
%
%         srange1 = t * (ScaleParameter(i) + PointCloud.ResolutionLocal(Keypoint.Location(KeypointsAtLevel(j))));
%         srange2 = t * (ScaleParameter(i+1) + PointCloud.ResolutionLocal(Keypoint.Location(KeypointsAtLevel(j))));
%
%         [Neigh1, ~] = rangesearch(KeypointTree, CurrentPoint, srange1);
%         [Neigh2, ~] = rangesearch(KeypointTree, CurrentPoint, srange2);
%
%         Neigh1 = Neigh1{:};
%         Neigh2 = Neigh2{:};
%
%
%         ZeroSphere1 = ismember(Neigh1, ZeropointAtLevel);
%         ZeroSphere2 = ismember(Neigh2, ZeropointAtLevel);
%
%         if ~nnz(ZeroSphere1) && nnz(ZeroSphere2)
%             NumNewKeypoint = NumNewKeypoint + 1;
%             NewKeypoint.Level(NumNewKeypoint,1) = i;
%             NewKeypoint.Location(NumNewKeypoint, 1) = Keypoint.Location(KeypointsAtLevel(j));
%             NewKeypoint.Scale(NumNewKeypoint, 1) = srange1 / t;
%             NewKeypoint.Count = NewKeypoint.Count + 1;
%         end
%
%     end
%
%     waitbar(i/NumLevels, WaitBar, sprintf('Checking Keypoint %i of %i', i, NumLevels));
%
% end






% Remove zero locations
ZeroLogic = NewKeypoint.Level == 0;
NewKeypoint.Level(ZeroLogic) = [];
NewKeypoint.Location(ZeroLogic) = [];
NewKeypoint.Scale(ZeroLogic) = [];
NewKeypoint.Sign(ZeroLogic) = [];
NewKeypoint.DoG(ZeroLogic) = [];





end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







































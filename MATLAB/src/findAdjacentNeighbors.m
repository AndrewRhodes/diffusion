% Andrew Rhodes
% ASEL 2017
% find the adjacent neighbors for each point
%
% input[PointCloudIn]: structure with parameters 'Location',
%                   'LocationCount', 'Face', 'FaceCount', 'Resolution'
%
% input[[Type]; string, 'connectivity', 'distance'
%
% output[Neighbors]: cell size[nx1]



function [Neighbors, NeighborFaces, PointCloudIn] = findAdjacentNeighbors(PointCloudIn)

if ~isfield(PointCloudIn, 'Location')
    erorr('PointCloudIn must contain field ''Location''')
end

if ~isfield(PointCloudIn, 'Face')
    erorr('PointCloudIn must contain field ''Face''')
end

if ~isfield(PointCloudIn, 'LocationCount')
    warning('Adding field ''LocationCount'' to PointCloudIn')
    PointCloudIn.LocationCount = length(PointCloudIn.Location);
end

if ~isfield(PointCloudIn, 'FaceCount')
    warning('Adding field ''FaceCount'' to PointCloudIn')
    PointCloudIn.FaceCount = length(PointCloudIn.Face);
end


Neighbors.Connect = cell(PointCloudIn.LocationCount,1);
NeighborFaces.Connect = cell(PointCloudIn.LocationCount,1);


for i = 1 : PointCloudIn.FaceCount
    Neighbors.Connect{PointCloudIn.Face(i,1)} = [Neighbors.Connect{PointCloudIn.Face(i,1)} [PointCloudIn.Face(i,2) PointCloudIn.Face(i,3)]];
    Neighbors.Connect{PointCloudIn.Face(i,2)} = [Neighbors.Connect{PointCloudIn.Face(i,2)} [PointCloudIn.Face(i,3) PointCloudIn.Face(i,1)]];
    Neighbors.Connect{PointCloudIn.Face(i,3)} = [Neighbors.Connect{PointCloudIn.Face(i,3)} [PointCloudIn.Face(i,1) PointCloudIn.Face(i,2)]];
    
    NeighborFaces.Connect{PointCloudIn.Face(i,1)} = [NeighborFaces.Connect{PointCloudIn.Face(i,1)}, i];
    NeighborFaces.Connect{PointCloudIn.Face(i,2)} = [NeighborFaces.Connect{PointCloudIn.Face(i,2)}, i];
    NeighborFaces.Connect{PointCloudIn.Face(i,3)} = [NeighborFaces.Connect{PointCloudIn.Face(i,3)}, i];
    
end

    
Neighbors.Connect = cellfun(@unique, Neighbors.Connect, 'UniformOutput', 0);
NeighborFaces.Connect = cellfun(@unique, NeighborFaces.Connect, 'UniformOutput', 0);


% Calculate local resolution of point cloud. And pass back out
PointCloudIn = findLocalResolution(PointCloudIn, Neighbors.Connect);


KDTree = KDTreeSearcher(PointCloudIn.Location, 'Distance', 'euclidean');

NeighborFaces.Distance = cell(PointCloudIn.LocationCount,1);
Neighbors.Distance = cell(PointCloudIn.LocationCount,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distance based neighbor search in loop format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 1 : PointCloudIn.LocationCount
%     
%     [Neigh, ~] = rangesearch(KDTree, PointCloudIn.Location(i,:), sqrt(3)*PointCloudIn.ResolutionLocal(i,1));
%     
%     Neigh{1}(1) = [];
%     
%     Neighbors.Distance{i} = [Neigh{1}];
%     
%     [~, IndexA, ~] = intersect(PointCloudIn.Face,[Neigh{1}]);
%     
%     NeighborFaces.Distance{i} = IndexA;
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distance based neighbor search in vectorized format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


searcher = @(Tree, Point, Radius) rangesearch(KDTree, Point, Radius);
remover = @(neighbor) (neighbor(2:end));

Neig = cellfun(searcher, repmat({KDTree}, PointCloudIn.LocationCount, 1),...
        mat2cell(PointCloudIn.Location, ones(PointCloudIn.LocationCount,1), 3),...
        num2cell(sqrt(3)*PointCloudIn.ResolutionLocal)); % sqrt(3) because in R^3. 

Neighbors.Distance = cellfun(remover, Neig, 'Uni', 0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Neighbor Faces based on Distance
% Takes a long time


% [~, NeighborFaces.Distance, ~] = cellfun(@intersect,...
%     repmat(mat2cell(PointCloudIn.Face, PointCloudIn.FaceCount, 3), PointCloudIn.LocationCount,1),...
%     Neighbors.Distance, 'Uni', 0);





end



















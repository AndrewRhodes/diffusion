% Andrew Rhodes
% ASEL 2017
% find the adjacent neighbors for each point
%
% input[faces]: size[mx1]
% input[vertices]: size[nx1]
% 
% output[Neighbors]: cell size[nx1]



function Neighbors = findAdjacentNeighbors(PointCloudIn)

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





Neighbors = cell(PointCloudIn.LocationCount,1);

for i = 1 : PointCloudIn.FaceCount
    Neighbors{PointCloudIn.Face(i,1)} = [Neighbors{PointCloudIn.Face(i,1)} [PointCloudIn.Face(i,2) PointCloudIn.Face(i,3)]];
    Neighbors{PointCloudIn.Face(i,2)} = [Neighbors{PointCloudIn.Face(i,2)} [PointCloudIn.Face(i,3) PointCloudIn.Face(i,1)]];
    Neighbors{PointCloudIn.Face(i,3)} = [Neighbors{PointCloudIn.Face(i,3)} [PointCloudIn.Face(i,1) PointCloudIn.Face(i,2)]];
end

Neighbors = cellfun(@unique, Neighbors, 'UniformOutput', 0);


end
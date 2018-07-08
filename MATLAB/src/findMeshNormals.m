% Andrew Rhodes
% ASEL
% November 2015
% -----
% input[PointCloudIn]: Structure that must at least contain 
%                       .Location & .Face
% 
% output[PointCloudIn]: Same as input structure with addition of .Normal
% -----
% The program calculates a normal for each mesh face. Then the vertex
% normals are a weighted sum of their connecting faces' normals.


function PointCloudIn = findMeshNormals(PointCloudIn)

Fields = fieldnames(PointCloudIn);

if ~strcmp(Fields, 'Location')
    error('Input Structure PointCloudIn must contain a field .Location')
elseif ~strcmp(Fields, 'Face')
    error('Input Structure PointCloudIn must contain a field .Face')
end

if size(PointCloudIn.Location, 2) ~= 3
    error('Input Structure PointCloudIn.Location must have size [mx3].')
end

if size(PointCloudIn.Face, 2) ~= 3
    error('Input Structure PointCloudIn.Face must have size [nx3].')
end

if ~strcmp(Fields, 'Count')
   PointCloudIn.Count = size(PointCloudIn.Location, 1); 
end



% Get all edge vectors
e1 = PointCloudIn.Location(PointCloudIn.Face(:,1),:) - PointCloudIn.Location(PointCloudIn.Face(:,2),:);
e2 = PointCloudIn.Location(PointCloudIn.Face(:,2),:) - PointCloudIn.Location(PointCloudIn.Face(:,3),:);
e3 = PointCloudIn.Location(PointCloudIn.Face(:,3),:) - PointCloudIn.Location(PointCloudIn.Face(:,1),:);

% Normalize edge vectors
e1 = bsxfun(@rdivide, e1, sqrt(sum(e1.^2,2)));
e2 = bsxfun(@rdivide, e2, sqrt(sum(e2.^2,2)));
e3 = bsxfun(@rdivide, e3, sqrt(sum(e3.^2,2)));

% Calculate Angle of face seen from vertices
Angle = [acos(dot(e1',-e3')); acos(dot(e2',-e1')); acos(dot(e3',-e2'))]';

% Calculate normal of face
Normal=cross(e1,e3);

% NumFaces = size(FacesIn,1);
% NumLocations = size(VerticesIn,1);
% NormalsOut = zeros(NumLocations, 3);
PointCloudIn.Normal = zeros(PointCloudIn.Count, 3);


% Calculate Vertice Normals
PointCloudIn.Normal(PointCloudIn.Face(:,1),:) = bsxfun(@plus, PointCloudIn.Normal(PointCloudIn.Face(:,1),:), bsxfun(@times, Normal, Angle(:,1)));
PointCloudIn.Normal(PointCloudIn.Face(:,2),:) = bsxfun(@plus, PointCloudIn.Normal(PointCloudIn.Face(:,2),:), bsxfun(@times, Normal, Angle(:,2)));
PointCloudIn.Normal(PointCloudIn.Face(:,3),:) = bsxfun(@plus, PointCloudIn.Normal(PointCloudIn.Face(:,3),:), bsxfun(@times, Normal, Angle(:,3)));

% Normals point towards the center
PointCloudIn.Normal = bsxfun(@rdivide, PointCloudIn.Normal, sqrt(sum(PointCloudIn.Normal.^2,2)));

% flip normals to point outward
PointCloudIn.Normal = - PointCloudIn.Normal;



% Clean up instances of Normals of NaN.
NanNormal = find(isnan(PointCloudIn.Normal(:,1)))

for i = 1 : length(NanNormal)
    
    [row, ~] = ind2sub([PointCloudIn.FaceCount, 3], find(PointCloudIn.Face == NanNormal(i)));
    
    NeighborPoints = unique(PointCloudIn.Face(row,:));
    NeighborPoints(NeighborPoints == NanNormal(i)) = [];
    NeighborNormals = PointCloudIn.Normal(NeighborPoints,:);
    NeighborNormals = reshape(NeighborNormals(~isnan(NeighborNormals)),[],3);
    PointCloudIn.Normal(NanNormal(i),:) = sum(NeighborNormals,1)/size(NeighborNormals,1);
    
    
end






end







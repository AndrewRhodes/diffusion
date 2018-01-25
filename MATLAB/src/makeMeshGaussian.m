% Andrew Rhodes
% WVU
% Jan. 2018
%
% makeMeshGaussian makes the approximate Gaussian on triangulated mesh
% surfaces
%
% input[PointCloud]: structure must have elements Location of size [mx3]
% and Face size [nx3]
% input[AverageEdgeLength]: average edge length of mesh
% input[SearchRange]: the width of the gaussian, usually 2*sqrt(2)*AverageEdgeLength
% input[NormalizeLogic]:
%
% output[]:



function G = makeMeshGaussian(PointCloud, AverageEdgeLength, SearchRange, NormalizeLogic)

if ~isfield(PointCloud, {'Location','Face', 'FaceArea'})
    error('PointCloud structure needs elements Location, Face, and FaceArea')
end

if ~isfield(PointCloud, 'FaceCount')
    PointCloud.FaceCount = length(PointCloud.Face);
end
if ~isfield(PointCloud, 'LocationCount')
    PointCloud.LocationCount = length(PointCloud.Location);
end



KDTree = KDTreeSearcher(PointCloud.Location);


if ~isfield(PointCloud, 'VertexArea')
    PointCloud.VertexArea = accumarray( reshape(PointCloud.Face, 3*PointCloud.FaceCount ,1), repmat(PointCloud.FaceArea,3,1), [PointCloud.LocationCount, 1] )/3;
end

% Initialize Vector Sizes
Position1 = zeros(25 * PointCloud.LocationCount,1);
Position2 = zeros(25 * PointCloud.LocationCount,1);
Weights = zeros(25 * PointCloud.LocationCount,1);



count = 1;
maxcount = 0;
WaitBar = waitbar(0, sprintf('Building Gaussian Vertex %i of %i', 0, PointCloud.LocationCount-1));

for i = 1 : PointCloud.LocationCount
    
    
    [Neigh, ~] = rangesearch(KDTree, PointCloud.Location(i,:), SearchRange);
    
    CurrentNeighbors = Neigh{:};
    
    NumNeighbors = length(CurrentNeighbors);
    
    maxcount = maxcount + NumNeighbors;
    
    Position1(count:maxcount) = i * ones(NumNeighbors,1);
    Position2(count:maxcount) = CurrentNeighbors;
    
    Weights(count:maxcount) = 1/(2*pi*AverageEdgeLength^2) * PointCloud.VertexArea(i) * PointCloud.VertexArea(CurrentNeighbors) .* exp(- sum(bsxfun(@minus,PointCloud.Location(i,:), PointCloud.Location(CurrentNeighbors,:)).^2,2) / (2*AverageEdgeLength^2)  );
    
    
    count = maxcount + 1;
    waitbar(i/PointCloud.LocationCount, WaitBar, sprintf('Building Gaussian Vertex %i of %i', i, PointCloud.LocationCount-1));
    
end

waitbar(i/PointCloud.LocationCount, WaitBar, sprintf('Building Complete'));
close(WaitBar)



Weights(Position1==0) = [];
Position2(Position1==0) = [];
Position1(Position1==0) = [];

GA = sparse(Position1, Position2, Weights, PointCloud.LocationCount, PointCloud.LocationCount);

B = sparse(1:PointCloud.LocationCount, 1:PointCloud.LocationCount, 1./PointCloud.VertexArea) * GA;


if NormalizeLogic
    G = B ./ max(sum(B,2));
else
    G = B;
end


end







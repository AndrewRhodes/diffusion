% Andrew Rhodes
% WVU
% Jan. 2018
%
% makeMeshGaussian makes the approximate Gaussian on triangulated mesh
% surfaces
% 
% The Gaussian std is assumed to be the AverageEdgeLength
% 
% With Matlab, it can only use Eucluidean distance between neighboring
% verticies
%
% input[PointCloud]: structure must have elements Location of size [nx3]
% and Face size [mx3]
% input[Sigma]: standard deviation for the gaussian
% input[NumSigma]: the support region for the gaussian
% input[NormalizeMethod]: Normalize the rows to unity
%
% output[G]: The discrete mesh Gaussian, sparse, size[nxn]


% Thoughts: Could AverageEdgeLength be any arbitrary scalar?
% Thoughts: Could make this adaptive by changing the AverageEdgeLength to
%           be a function of the edge length of the first-ring, or 
%           second-ring neighbors.



function G = makeExplicitGaussian(PointCloud, Sigma, NumSigma)

% Determine if 
if ~isfield(PointCloud, {'Location','Face', 'FaceArea'})
    error('PointCloud structure needs elements Location, Face, and FaceArea')
end

if ~isfield(PointCloud, 'FaceCount')
    PointCloud.FaceCount = length(PointCloud.Face);
end
if ~isfield(PointCloud, 'LocationCount')
    PointCloud.LocationCount = length(PointCloud.Location);
end


% Find the vertex area weight if not already calculated.
if ~isfield(PointCloud, 'VertexArea')
    warning('Calculating PointCloud.VertexArea')
    PointCloud.VertexArea = accumarray( reshape(PointCloud.Face, 3*PointCloud.FaceCount ,1), repmat(PointCloud.FaceArea,3,1), [PointCloud.LocationCount, 1] )/3;
end



% Build a KDTree for searching nearby neighbors.
KDTree = KDTreeSearcher(PointCloud.Location, 'distance', 'euclidean');



% Find all neighbor points within a radius
[Neigh, ~] = rangesearch(KDTree, PointCloud.Location, Sigma*NumSigma);

% Calculate the number of neighbors per each vertex
NumNeighbors = cellfun(@length, Neigh);

% Create indexing vectors for sparse matrix below
Position2 = horzcat(Neigh{:})';

MakePos1 = @(Ind, NumNei) Ind*ones(NumNei,1);
Pos1 = cellfun(MakePos1, num2cell(1:PointCloud.LocationCount)', num2cell(NumNeighbors), 'uniformoutput',0);
Position1 = vertcat(Pos1{:});


% Place critical items in cells for use of cellfun
VerNeighArea = cell(PointCloud.LocationCount,1);
PCLocNeighCell = cell(PointCloud.LocationCount,1);
for i = 1 : PointCloud.LocationCount
%    Concatenate area of neighbor vertices in cells
   VerNeighArea{i} = PointCloud.VertexArea(Neigh{i});
   
%    Concatenate location of neighbor vertices in cells
   PCLocNeighCell{i} = PointCloud.Location(Neigh{i},:);
end


% The ratio of the center vertex area with it's neighbor vertex areas
AreaRatio = cellfun(@rdivide, VerNeighArea, num2cell(PointCloud.VertexArea), 'uniformoutput',0); 
% AreaRatio = cellfun(@rdivide, VerNeighArea, VerNeighArea, 'uniformoutput',0);

% Place PointCloud.Location into a cell for later use in cellfun
PCLocCell = mat2cell(PointCloud.Location,ones(PointCloud.LocationCount,1),3);


% Define the function for the Gaussian exponential part for use in cellfun
ExpDiffFunc = @(CentVert, NeighVert, AreaRatio, sigma) exp( - AreaRatio .* sqrt( sum( bsxfun( @minus, CentVert, NeighVert).^2, 2) ) ./ (2*sigma^2) );
% Find the exponential part of the Gaussian using cellfun
WeightExp = cellfun(ExpDiffFunc, PCLocCell, PCLocNeighCell, AreaRatio, num2cell(Sigma*ones(PointCloud.LocationCount,1)),'uniformoutput',0);


% Multiple the Gaussian constant with the exponential
% ConstExpWeight = cellfun(@times, num2cell(1/((2*pi)^(3/2)*Sigma^3)*ones(PointCloud.LocationCount,1)), AreaExpWeight, 'uniformoutput',0);

% Reshape into column vector for use in sparse matrix constuction
Weights = vertcat(WeightExp{:});

% Weights = vertcat(AreaExpWeight{:});

GA = sparse(Position1, Position2, Weights, PointCloud.LocationCount, PointCloud.LocationCount);



G = bsxfun(@rdivide, GA , sum(GA,2));



end








%% Second incarnation. This is with the weight infront of the Gaussian
% function G = makeExplicitGaussian(PointCloud, AverageEdgeLength, SearchRange, NormalizeMethod)
% 
% % Determine if 
% if ~isfield(PointCloud, {'Location','Face', 'FaceArea'})
%     error('PointCloud structure needs elements Location, Face, and FaceArea')
% end
% 
% if ~isfield(PointCloud, 'FaceCount')
%     PointCloud.FaceCount = length(PointCloud.Face);
% end
% if ~isfield(PointCloud, 'LocationCount')
%     PointCloud.LocationCount = length(PointCloud.Location);
% end
% 
% % Build a KDTree for searching nearby neighbors.
% KDTree = KDTreeSearcher(PointCloud.Location);
% 
% % Find the vertex area weight if not already calculated.
% if ~isfield(PointCloud, 'VertexArea')
%     warning('Calculating PointCloud.VertexArea')
%     PointCloud.VertexArea = accumarray( reshape(PointCloud.Face, 3*PointCloud.FaceCount ,1), repmat(PointCloud.FaceArea,3,1), [PointCloud.LocationCount, 1] )/3;
% end
% 
% 
% % Find all neighbor points within a radius
% [Neigh, ~] = rangesearch(KDTree, PointCloud.Location, SearchRange);
% 
% % Calculate the number of neighbors per each vertex
% NumNeighbors = cellfun(@length, Neigh);
% 
% % Create indexing vectors for sparse matrix below
% Position2 = horzcat(Neigh{:})';
% 
% MakePos1 = @(Ind, NumNei) Ind*ones(NumNei,1);
% Pos1 = cellfun(MakePos1, num2cell(1:PointCloud.LocationCount)', num2cell(NumNeighbors), 'uniformoutput',0);
% Position1 = vertcat(Pos1{:});
% 
% 
% % Place critical items in cells for use of cellfun
% VerNeighArea = cell(PointCloud.LocationCount,1);
% PCLocNeighCell = cell(PointCloud.LocationCount,1);
% for i = 1 : PointCloud.LocationCount
% %    Concatenate area of neighbor vertices in cells
%    VerNeighArea{i} = PointCloud.VertexArea(Neigh{i});
%    
% %    Concatenate location of neighbor vertices in cells
%    PCLocNeighCell{i} = PointCloud.Location(Neigh{i},:);
% end
% 
% 
% % The product of the center vertex area with it's neighbor vertex areas
% AreaProd = cellfun(@times, num2cell(PointCloud.VertexArea), VerNeighArea, 'uniformoutput',0);
% % AreaProd = VerNeighArea;
% 
% % Place PointCloud.Location into a cell for later use in cellfun
% PCLocCell = mat2cell(PointCloud.Location,ones(PointCloud.LocationCount,1),3);
% 
% % Define the function for the Gaussian exponential part for use in cellfun
% ExpDiffFunc = @(CentVert, NeighVert, sigma) exp( - sqrt( sum( bsxfun( @minus, CentVert, NeighVert).^2, 2) ) ./ (2*sigma^2) );
% % Find the exponential part of the Gaussian using cellfun
% WeightExp = cellfun(ExpDiffFunc, PCLocCell, PCLocNeighCell, num2cell(AverageEdgeLength*ones(PointCloud.LocationCount,1)),'uniformoutput',0);
% 
% % Mulitple the area weight part with the exponetial of the Gaussian
% AreaExpWeight = cellfun(@times, AreaProd, WeightExp, 'uniformoutput',0);
% % AreaExpWeight = WeightExp;
% 
% 
% % Multiple the Gaussian constant with the other terms
% ConstAreaExpWeight = cellfun(@times, num2cell(1/(2*pi*AverageEdgeLength^2)*ones(PointCloud.LocationCount,1)), AreaExpWeight, 'uniformoutput',0);
% % ConstAreaExpWeight = cellfun(@times, num2cell(ones(PointCloud.LocationCount,1)), AreaExpWeight, 'uniformoutput',0);
% 
% % Reshape into column vector for use in sparse matrix constuction
% Weights = vertcat(ConstAreaExpWeight{:});
% 
% % Weights = vertcat(AreaExpWeight{:});
% 
% GA = sparse(Position1, Position2, Weights, PointCloud.LocationCount, PointCloud.LocationCount);
% 
% B = sparse(1:PointCloud.LocationCount, 1:PointCloud.LocationCount, 1./PointCloud.VertexArea) * GA;
% 
% 
% if NormalizeMethod == 1
%     G = bsxfun(@rdivide, GA , sum(GA,2));
% elseif NormalizeMethod == 2
%     G = bsxfun(@rdivide, B , sum(B,2));
%     G = bsxfun(@rdivide, G , sum(G,2));
% end
% 
% 
% end




%% Original code for constructing the Gaussian. 
% For loop is too slow.

% Weights = zeros(length(Position2),1);
% 
% count = 1;
% maxcount = 0;
% WaitBar = waitbar(0, sprintf('Building Gaussian Vertex %i of %i', 0, PointCloud.LocationCount-1));
% 
% for i = 1 : PointCloud.LocationCount
%          
%     maxcount = maxcount + NumNeighbors(i);
%     
%     Weights(count:maxcount) = 1/(2*pi*AverageEdgeLength^2) * PointCloud.VertexArea(i) * PointCloud.VertexArea(Neigh{i,1}) .* exp(- sum(bsxfun(@minus,PointCloud.Location(i,:), PointCloud.Location(Neigh{i,1},:)).^2,2) / (2*AverageEdgeLength^2)  );
%     
%     count = maxcount + 1;
%     waitbar(i/PointCloud.LocationCount, WaitBar, sprintf('Building Gaussian Vertex %i of %i', i, PointCloud.LocationCount-1));
%     
% end
% 
% waitbar(i/PointCloud.LocationCount, WaitBar, sprintf('Building Complete'));
% close(WaitBar)
% 
% Weights(Position1==0) = [];
% Position2(Position1==0) = [];
% Position1(Position1==0) = [];





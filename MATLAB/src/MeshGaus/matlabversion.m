


global ProjectRoot;
FileLocationModel = strcat(ProjectRoot,'/models/object/');
ModelFolder = 'armadillo/';
Model = 'Armadillo_e1_100000';
FileNameModelOff = strcat(ModelFolder,Model,'.off');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');

Sigma = h
NumSigma = opt.rho

[PointCloud.Location, PointCloud.Face, PointCloud.Normal, PointCloud.Signal]...
                = read_ply_all_elements( fullfile( FileLocationModel, FileNameModelPly ) );

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
PointCloud = findMeshResolution(PointCloud, 'Model');



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
PCLocNeighCell = cell(PointCloud.LocationCount,1);
for i = 1 : PointCloud.LocationCount
%    Concatenate location of neighbor vertices in cells
   PCLocNeighCell{i} = PointCloud.Location(Neigh{i},:);
end


% The ratio of the center vertex area with it's neighbor vertex areas
% AreaRatio = cellfun(@rdivide, VerNeighArea, VerNeighArea, 'uniformoutput',0);

% Place PointCloud.Location into a cell for later use in cellfun
PCLocCell = mat2cell(PointCloud.Location,ones(PointCloud.LocationCount,1),3);


% Define the function for the Gaussian exponential part for use in cellfun
ExpDiffFunc = @(CentVert, NeighVert, sigma) (4/(pi*sigma^4)) * exp( - sum( bsxfun( @minus, CentVert, NeighVert).^2, 2) ./ (sigma^2) );
% Find the exponential part of the Gaussian using cellfun
WeightExp = cellfun(ExpDiffFunc, PCLocCell, PCLocNeighCell, num2cell(Sigma*ones(PointCloud.LocationCount,1)),'uniformoutput',0);


% Reshape into column vector for use in sparse matrix constuction
Weights = vertcat(WeightExp{:});

% Weights = vertcat(AreaExpWeight{:});

GA = sparse(Position1, Position2, Weights, PointCloud.LocationCount, PointCloud.LocationCount);


Gm = bsxfun(@rdivide, GA , sum(GA,2));


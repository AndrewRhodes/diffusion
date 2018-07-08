% Andrew Rhodes
% WVU
% Jan. 2018
%
% findKeyPoints finds the keypoint locations and their scale on 2D images
% or 3D meshes. 
%
% input[DoG]: Difference of Gaussian of signals over which to search
% input[ScaleParameter]: list of scale parameters accompanying the DoG
% input[Neighbors]: comparison neighbors for determining local extrema
% input[NumPoints]: the number of points that are being searched
% input[varargin]: dimension of data, either 2 or 3. Default is 3
%
% output[KeyPoint]: structure with elements .Location and .Scale






function Keypoint = findKeypoint2(DoG, ScaleParameter, Neighbors, varargin)


if nargin < 4
    Dimension = 3;
elseif nargin == 4
    Dimension = varargin{1};
elseif nargin > 4
    error('A maximum of 5 arguments is permitted.')
end


if ~iscell(Neighbors)
    error('Neighbors must be in cell format')
end

% % % Pre-computation checks and setup
[NumVertices, NumDoGLevels] = size(DoG);

% Number of Neighbors per vertex
NeighborLength = cell2mat(cellfun(@length, Neighbors, 'Uni', 0));


% Number of Above/Below Neighbors per vertex. Equivalent Statements
% NeighborAboveBelowLength = cell2mat(cellfun(@length, NeighborAboveBelowCell, 'Uni', 0));
NeighborAboveBelowLength = NeighborLength + 1;
sumNeighborAboveBelowLength = sum(NeighborAboveBelowLength);

CurrentLevelDoGNeighbors = zeros(sum(NeighborLength), NumDoGLevels);
AboveLevelDoGNeighbors = zeros(sumNeighborAboveBelowLength, NumDoGLevels);
BelowLevelDoGNeighbors = zeros(sumNeighborAboveBelowLength, NumDoGLevels);
clear sumNeighborAboveBelowLength

NumVerticesCell = num2cell(1:NumVertices)';
% Concatenate Above/Below Neighbor Indices into array -include center vertex
NeighborAboveBelowMat = cell2mat(cellfun(@(N,N2) [N';N2], Neighbors, NumVerticesCell, 'Uni',0));

% Convert Neighbors into 1d array
NeighborMat = cell2mat(cellfun(@(N) N',Neighbors, 'Uni', 0));

CurrentLevelDoGNeighbors(:,1) = DoG(NeighborMat,1);
CurrentLevelDoGNeighbors(:,NumDoGLevels) = DoG(NeighborMat,NumDoGLevels);
for j = 2 : NumDoGLevels - 1
    CurrentLevelDoGNeighbors(:,j) = DoG(NeighborMat,j);
    AboveLevelDoGNeighbors(:,j) = DoG(NeighborAboveBelowMat,j+1);
    BelowLevelDoGNeighbors(:,j) = DoG(NeighborAboveBelowMat,j-1);
end
clear NeighborMat NeighborAboveBelowMat j

% Convert Neighbor DoG to cell.
CurrentLevelDoGNeighborsCell = mat2cell(CurrentLevelDoGNeighbors, NeighborLength, ones(NumDoGLevels,1));
clear CurrentLevelDoGNeighbors NeighborLength

% Convert Above/Below Neighbor DoG to cell. Reassign first/last columns.
AboveLevelDoGNeighborsCell = mat2cell(AboveLevelDoGNeighbors, NeighborAboveBelowLength, ones(NumDoGLevels,1));
clear AboveLevelDoGNeighbors

AboveLevelDoGNeighborsCell(:,1) = NumVerticesCell;
AboveLevelDoGNeighborsCell(:,NumDoGLevels) = NumVerticesCell;

BelowLevelDoGNeighborsCell = mat2cell(BelowLevelDoGNeighbors, NeighborAboveBelowLength, ones(NumDoGLevels,1));
clear BelowLevelDoGNeighbors NeighborAboveBelowLength

BelowLevelDoGNeighborsCell(:,1) = NumVerticesCell;
BelowLevelDoGNeighborsCell(:,NumDoGLevels) = NumVerticesCell;
clear NumVerticesCell

% Define Anon Function to test >/< DoG values
FindMax = @(DoG, CurrentLevelDoG, AboveLevelDoG, BelowLevelDog) ...
(all(DoG > CurrentLevelDoG) & all(DoG > AboveLevelDoG) & all(DoG > BelowLevelDog)) | ...
(all(DoG < CurrentLevelDoG) & all(DoG < AboveLevelDoG) & all(DoG < BelowLevelDog));

KeypointIndices = cellfun(FindMax, num2cell(DoG), CurrentLevelDoGNeighborsCell, AboveLevelDoGNeighborsCell, BelowLevelDoGNeighborsCell, 'Uni', 0);
clear CurrentLevelDoGNeighborsCell AboveLevelDoGNeighborsCell BelowLevelDoGNeighborsCell FindMax DoG
 
[Location, Level] = find(cell2mat(KeypointIndices));

% Remove first and last column of DoG Levels
Logic = (Level ~= 1 & Level ~= NumDoGLevels);

% Create Keypoint Structure
Keypoint.Location = Location(Logic);
Keypoint.Level = Level(Logic);
Keypoint.Scale = ScaleParameter(Keypoint.Level);
Keypoint.Count = length(Keypoint.Location);




end




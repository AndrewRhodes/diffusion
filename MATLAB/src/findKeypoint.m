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





function Keypoint = findKeypoint(DoG, PointCloud, ScaleParameter, Neighbors, Method, CompareMethod)

%
% if nargin < 4
%     Dimension = 3;
% elseif nargin == 4
%     Dimension = varargin{1};
% elseif nargin > 4
%     error('A maximum of 5 arguments is permitted.')
% end

% % % Pre-computation checks and setup
isSign = 0;

if strcmpi(Method, 'Old')
    
    [NumVertices, NumDoGLevels] = size(DoG);
    
    if ~iscell(Neighbors)
        error('Neighbors must be in cell format')
    end
    
    
    MaxNumberFeatures = round(NumDoGLevels * NumVertices / 100);
    
    
    % Initialize space for the keypoints
    
    
    Keypoint.Scale = zeros(MaxNumberFeatures,1);
    Keypoint.Location = zeros(MaxNumberFeatures,3);
    Keypoint.LocationIndex = zeros(MaxNumberFeatures,1);
    Keypoint.Normal = zeros(MaxNumberFeatures,3);
    Keypoint.Level = zeros(MaxNumberFeatures,1);
    Keypoint.Sign = zeros(MaxNumberFeatures,1);
    
    NumFeatures = 0;
    
    isFeaturePoint = 0;
    
    WaitBar = waitbar(0, sprintf('Checking Vertex %i of %i', 0, NumVertices));
    
    for i = 1 : NumVertices
        
        
        %     if Dimension == 2
        %         CurrentNeighbors = sub2ind([ImageSize(1),ImageSize(2)], CurrentNeighbors(:,1), CurrentNeighbors(:,2));
        %     elseif Dimension == 3
        CurrentNeighbors = Neighbors{i,1};
        %     else
        %         error('Dimension must be 2 or 3.')
        %     end
        
        
        for j = 2 : NumDoGLevels - 2
            
            CurrentValue = DoG(i,j);
            
            SurroundingValuesBelow = DoG([i,CurrentNeighbors],j-1);
            SurroundingValuesCurrent = DoG(CurrentNeighbors,j);
            SurroundingValuesAbove = DoG([i,CurrentNeighbors],j+1);
            
            if strcmpi(CompareMethod, '>')
                
                if all(CurrentValue > SurroundingValuesCurrent)
                    if all(CurrentValue > SurroundingValuesBelow)
                        if all(CurrentValue > SurroundingValuesAbove)
                            isFeaturePoint = 1;
                        end
                    end
                end
                
            elseif strcmpi(CompareMethod, '<>') || strcmpi(CompareMethod, '><')
                
                if all(CurrentValue > SurroundingValuesCurrent)
                    if all(CurrentValue > SurroundingValuesBelow)
                        if all(CurrentValue > SurroundingValuesAbove)
                            isFeaturePoint = 1;
                            isSign = 1;
                        end
                    end
                elseif all(CurrentValue < SurroundingValuesCurrent)
                    if all(CurrentValue < SurroundingValuesBelow)
                        if all(CurrentValue < SurroundingValuesAbove)
                            isFeaturePoint = 1;
                            isSign = -1;
                        end
                    end
                end
                
            elseif strcmpi(CompareMethod, '<')
                
                if all(CurrentValue < SurroundingValuesCurrent)
                    if all(CurrentValue < SurroundingValuesBelow)
                        if all(CurrentValue < SurroundingValuesAbove)
                            isFeaturePoint = 1;
                        end
                    end
                end
                
            elseif strcmpi(CompareMethod, 'special')
                
                if CurrentValue < DoG(i, j-1)
                    if CurrentValue < DoG(i, j+1)
                        if all(CurrentValue < SurroundingValuesCurrent)
                            isFeaturePoint = 1;
                            isSign = -1;
                        end
                    end
                elseif CurrentValue > DoG(i, j-1)
                    if CurrentValue > DoG(i, j+1)
                        if all(CurrentValue > SurroundingValuesCurrent)
                            isFeaturePoint = 1;
                            isSign = 1;
                        end
                    end
                end
                
                
            else 
                error('findKeypoint:: incorrecte input for CompareMethod')
            end
            
            
            if isFeaturePoint
                NumFeatures = NumFeatures + 1;
                Keypoint.Scale(NumFeatures, 1) = ScaleParameter(j);
%                 Keypoint.Scale(NumFeatures, 1) = ScaleParameter(j) + PointCloud.ResolutionLocal(i,1);
                Keypoint.Location(NumFeatures, 1:3) = PointCloud.Location(i,:);
                Keypoint.Normal(NumFeatures, 1:3) = PointCloud.Normal(i,:);
                Keypoint.LocationIndex(NumFeatures, 1) = i;
                Keypoint.Level(NumFeatures, 1) = j;
                isFeaturePoint = 0;
                Keypoint.Sign(NumFeatures,1) = isSign;
                isSign = 0;
            end
            
            
        end
        waitbar(i/NumVertices, WaitBar, sprintf('Checking Vertex %i of %i', i, NumVertices));
        
    end
    
    waitbar(i/NumVertices, WaitBar, sprintf('Checking Vertex for Keypoint Complete'));
    close(WaitBar)
    
    
    % Remove empty elements from initilization
    ZeroLogic = (Keypoint.LocationIndex == 0);
    
    Keypoint.Scale(ZeroLogic,:) = [];
    Keypoint.Location(ZeroLogic,:) = [];
    Keypoint.LocationIndex(ZeroLogic,:) = [];
    Keypoint.Level(ZeroLogic,:) = [];
    Keypoint.Normal(ZeroLogic,:) = [];
    Keypoint.Sign(ZeroLogic,:) = [];
    Keypoint.Count = length(Keypoint.LocationIndex);
    
    
elseif strcmpi(Method, 'New')
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
    clear NeighborMat NeighborAboveBelowMat
    
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
    clear CurrentLevelDoGNeighborsCell AboveLevelDoGNeighborsCell BelowLevelDoGNeighborsCell FindMax
    
    [Location, Level] = find(cell2mat(KeypointIndices));
    
    % Remove first and last column of DoG Levels
    Logic = (Level ~= 1 & Level ~= NumDoGLevels);
    
    % Create Keypoint Structure
    Keypoint.LocationIndex = Location(Logic);
    Keypoint.Location = PointCloud.Location(Keypoint.LocationIndex,:);
    Keypoint.Normal = PointCloud.Normal(Keypoint.LocationIndex,:);
    Keypoint.Level = Level(Logic);
    Keypoint.Scale = ScaleParameter(Keypoint.Level);
    Keypoint.Count = length(Keypoint.LocationIndex);
    
else
    error('findKeypoint:: incorrecte input for Method')
end


end




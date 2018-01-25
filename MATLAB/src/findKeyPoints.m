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
% input[Dimension]: dimension of data, either 2 or 3 
%
% output[KeyPoint]: structure with elements .Location and .Scale



function KeyPoint = findKeyPoints(DoG, ScaleParameter, Neighbors, NumPoints, Dimension)


% % % Pre-computation checks and setup

if ~iscell(Neighbors)
    error('Neighbors must be in cell format')
end

if any(ismember(length(ScaleParameter), size(DoG)))
    NumLevels = length(ScaleParameter);
else    
    error('ScaleParameter and DoG must have the same number of elements')
end

if nargin < 4 
    error('Need at least 4 inputs')
elseif nargin < 5
    Dimension = 3;
end

if Dimension ~= 3 || Dimension ~= 2
    error('Dimension must be equal to 2 or 3')
end



MaxNumberFeatures = NumLevels * NumPoints / 100;


% Initialize space for the keypoints

KeyPoint.Scale = zeros(MaxNumberFeatures,1);
KeyPoint.Location = zeros(MaxNumberFeatures,1);


for i = 1 : NumPoints
    
    CurrentNeighbors = Neighbors{i,1};
    
    if Dimension == 2
        CurrentNeighbors = sub2ind([ImageSize(1),ImageSize(2)], CurrentNeighbors(:,1), CurrentNeighbors(:,2));
    end
    
    for j = 2 : NumLevels - 2
        
        CurrentValue = DoG(i,j);
        
        
        SurroundingValues_under = DoG([i;CurrentNeighbors],j-1);
        SurroundingValues_same = DoG(CurrentNeighbors,j);
        SurroundingValues_above = DoG([i;CurrentNeighbors],j+1);
        
        
        if all(CurrentValue > SurroundingValues_same)
            
            if all(CurrentValue > SurroundingValues_under)
                if all(CurrentValue > SurroundingValues_above)
                    NumFeatures = NumFeatures + 1;
                    % Maximum
                    KeyPoint.Scale(NumFeatures, 1) = ScaleParameter(j);
                    KeyPoint.Location(NumFeatures, 1) = i;
                end
            end
            
        elseif all(CurrentValue < SurroundingValues_same)
            
            if all(CurrentValue < SurroundingValues_under)
                if all(CurrentValue < SurroundingValues_above)
                    % Minimum
                    NumFeatures = NumFeatures + 1;
                    KeyPoint.Scale(NumFeatures, 1) = ScaleParameter(j);
                    KeyPoint.Location(NumFeatures, 1) = i;
                end
            end
        end
        
    end
    
end

KeyPoint.Scale(KeyPoint.Location == 0,:) = [];
KeyPoint.Location(KeyPoint.Location == 0) = [];




end




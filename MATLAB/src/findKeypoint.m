




function Keypoint = findKeypoint(DoG, ScaleParameter, Neighbors)

[NumVertices, NumDoGLevels] = size(DoG);


























NumFeatures = 0;

isFeaturePoint = 0;


WaitBar = waitbar(0, sprintf('Checking Vertex %i of %i', 0, NumVertices));

for i = 1 : NumVertices
    
    CurrentNeighbors = Neighbors{i,1};
    
    for j = 2 : NumDoGLevels - 2
        
        CurrentValue = DoG(i,j);
        
        SurroundingValuesBelow = DoG([i,CurrentNeighbors],j-1);
        SurroundingValuesCurrent = DoG(CurrentNeighbors,j);
        SurroundingValuesAbove = DoG([i,CurrentNeighbors],j+1);
        
        if all(CurrentValue > SurroundingValuesCurrent)
            if all(CurrentValue > SurroundingValuesBelow)
                if all(CurrentValue > SurroundingValuesAbove)
                    isFeaturePoint = 1;
                end
            end
        elseif all(CurrentValue < SurroundingValuesCurrent)
            if all(CurrentValue < SurroundingValuesBelow)
                if all(CurrentValue < SurroundingValuesAbove)
                    isFeaturePoint = 1;
                end
            end
        end
        
        if isFeaturePoint
            NumFeatures = NumFeatures + 1;
            Keypoint.Scale(NumFeatures, 1) = ScaleParameter(j);
            Keypoint.Location(NumFeatures, 1) = i;
            Keypoint.Level(NumFeatures, 1) = j;
            NumFeatures = NumFeatures + 1;
            isFeaturePoint = 0;
        end
        
        
        
%         SurroundingValues = [DoG([i,CurrentNeighbors],1,j-1);
%             DoG(CurrentNeighbors,1,j);
%             DoG([i,CurrentNeighbors],1,j+1)];


%         % Check Both if it is maximum or minimum
%         if all(CurrentValue > SurroundingValues)
%             NumFeatures = NumFeatures + 1;
%             % Maximum
%             FeaturePoint.Scale(NumFeatures, 1) = ScaleParamterSpatial(j);
%             FeaturePoint.Location(NumFeatures, 1) = i;
%         end
%         if all(CurrentValue < SurroundingValues)
%             % Minimum
%              NumFeatures = NumFeatures + 1;
%             FeaturePoint.Scale(NumFeatures, 1) = ScaleParamterSpatial(j);
%             FeaturePoint.Location(NumFeatures, 1) = i;
%         end
        
    end
    waitbar(i/NumVertices, WaitBar, sprintf('Checking Vertex %i of %i', i, NumVertices));

end

waitbar(i/NumVertices, WaitBar, sprintf('Checking Vertex for Keypoint Complete'));
close(WaitBar)


end
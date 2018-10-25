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





function [Keypoint, Zeropoint] = findKeypointNew(PointCloud, ScaleParameter, Neighbors, DoG, AbsDoG)


% % % Pre-computation checks and setup


[NumVertices, NumDoGLevels] = size(DoG);

if ~iscell(Neighbors)
    error('Neighbors must be in cell format')
end


MaxNumberFeatures = round(NumDoGLevels * PointCloud.LocationCount / 100);


% Initialize space for the keypoints
Keypoint.Scale = zeros(MaxNumberFeatures,1);
Keypoint.Location = zeros(MaxNumberFeatures,1);
Keypoint.Level = zeros(MaxNumberFeatures,1);
Keypoint.Sign = zeros(MaxNumberFeatures,1);
Keypoint.DoG = zeros(MaxNumberFeatures,1);
Keypoint.Count = 0;

Zeropoint.Scale = zeros(MaxNumberFeatures*100,1);
Zeropoint.Location = zeros(MaxNumberFeatures*100,1);
Zeropoint.Level = zeros(MaxNumberFeatures*100,1);
Zeropoint.Count = 0;

NumZeropoint = 0;
NumKeypoint = 0;
isSign = 0;
isKeypoint = 0;
isZeropoint = 0;

WaitBar = waitbar(0, sprintf('Checking Vertex %i of %i', 0, PointCloud.LocationCount));


for i = 1 : PointCloud.LocationCount
    
    CurrentNeighbors = Neighbors{i,1};

    for j = 2 : NumDoGLevels - 2
        
        CurrentDoG = DoG(i,j);
        CurrentAbsDoG = AbsDoG(i,j);
        
        SurroundingDoGCurrent = DoG(CurrentNeighbors,j);
        SurroundingAbsDoGCurrent = AbsDoG(CurrentNeighbors,j);
        
        
        if all(CurrentDoG > SurroundingDoGCurrent)
            isKeypoint = 1;
            isSign = 1;
        elseif all(CurrentDoG < SurroundingDoGCurrent)
            isKeypoint = 1;
            isSign = -1;
        end
        
        
%         if nnz( movprod(sign([CurrentDoG,SurroundingDoGCurrent']),2,2) < 0 ) > 4
%             isZeropoint = 1;
%         end
        
        if all(CurrentAbsDoG < SurroundingAbsDoGCurrent)
            isZeropoint = 1;
        end
            
        
        if isKeypoint
            NumKeypoint = NumKeypoint + 1;
            Keypoint.Scale(NumKeypoint, 1) = ScaleParameter(j) + PointCloud.ResolutionLocal(i,1);
            Keypoint.Location(NumKeypoint, 1) = i;
            Keypoint.Level(NumKeypoint, 1) = j;
            Keypoint.DoG(NumKeypoint, 1) = CurrentDoG;
            Keypoint.Sign(NumKeypoint, 1) = isSign;
            Keypoint.Count = Keypoint.Count + 1;
            
            isSign = 0;
            isKeypoint = 0;
        end
        
        if isZeropoint
            NumZeropoint = NumZeropoint + 1; 
            Zeropoint.Scale(NumZeropoint, 1) = ScaleParameter(j) + PointCloud.ResolutionLocal(i,1);
            Zeropoint.Location(NumZeropoint, 1) = i;
            Zeropoint.Level(NumZeropoint, 1) = j;
            Zeropoint.Count = Zeropoint.Count + 1;
            isZeropoint = 0;
        end
            
    end
    
    waitbar(i/NumVertices, WaitBar, sprintf('Checking Vertex %i of %i', i, NumVertices));
    
end

waitbar(i/NumVertices, WaitBar, sprintf('Checking Vertex for Keypoint Complete'));
close(WaitBar)


% Remove empty elements from initilization
ZeroLogic = (Keypoint.Location == 0);

Keypoint.Scale(ZeroLogic,:) = [];
Keypoint.Location(ZeroLogic,:) = [];
Keypoint.Level(ZeroLogic,:) = [];
Keypoint.Sign(ZeroLogic,:) = [];
Keypoint.DoG(ZeroLogic,:) = [];


ZeroLogic = (Zeropoint.Location == 0);

Zeropoint.Scale(ZeroLogic,:) = [];
Zeropoint.Location(ZeroLogic,:) = [];
Zeropoint.Level(ZeroLogic,:) = [];


end








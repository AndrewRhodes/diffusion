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





function DoGMaximum = findMHEdge(DoG, ScaleParameter, PointCloud, Neighbors, t_quant)

%
% if nargin < 4
%     Dimension = 3;
% elseif nargin == 4
%     Dimension = varargin{1};
% elseif nargin > 4
%     error('A maximum of 5 arguments is permitted.')
% end

% % % Pre-computation checks and setup



[NumVertices, NumDoGLevels] = size(DoG);

if ~iscell(Neighbors)
    error('Neighbors must be in cell format')
end


MaxNumberFeatures = round(NumDoGLevels * NumVertices / 100);


% Initialize space for the keypoints


% Keypoint.Scale = zeros(MaxNumberFeatures,1);
% Keypoint.Location = zeros(MaxNumberFeatures,1);
% Keypoint.Level = zeros(MaxNumberFeatures,1);
% Keypoint.LocationCell = cell(NumDoGLevels,1);

% DoGMaximum.Scale = zeros(MaxNumberFeatures*100,1);
% DoGMaximum.Location = zeros(MaxNumberFeatures*100,1);
% DoGMaximum.Level = zeros(MaxNumberFeatures*100,1);
% DoGMaximum.LocationCell = cell(NumDoGLevels,1);


% NumFeatures = 0;
% ep  = 0;
% DoGcell = cell(NumDoGLevels,1);

% isFeaturePoint = 0;
% isDoGMax = 0;

% Quants = quantile(DoG,[0.25,0.5,0.75]);

% WaitBar = waitbar(0, sprintf('Checking DoG %i of %i', 0, NumDoGLevels));

DoGSign = sign(DoG);

DoGSignProd = movprod(DoGSign,2,2);

[Location, Level] = find(DoGSignProd<0);

DoGMaximum.Location = Location;
DoGMaximum.Scale = ScaleParameter(Level) + PointCloud.ResolutionLocal(Location,1);
DoGMaximum.Level = Level;
% 
% 
% for j = 2 : NumDoGLevels - 2
%     
%     CurrentDoGLevel = DoG(:,j);
%     DoGOutOfBounds = (CurrentDoGLevel > Quants(3,j) + 1.5 * (Quants(3,j)-Quants(1,j))) | (CurrentDoGLevel < Quants(3,j) - 1.5 * (Quants(3,j)-Quants(1,j)));
%     CurrentDoGLevel(DoGOutOfBounds) = [];
%     
%     DoGcell{j,1}= CurrentDoGLevel;
% end
% 
% DoGcellSign = cellfun(@sign, DoGcell, 'Uni', 0);
%     
%     
%     CurrentQuant = quantile(CurrentDoGLevel, t_quant);
%     
%     DoGLessThanQuant = CurrentDoGLevel < CurrentQuant;
%     
%     Location = find(CurrentDoGLevel(DoGLessThanQuant));
%     
%     sp = ep + 1;
%     ep = sp + length(Location) - 1;
%     
%     DoGMaximum.Location(sp:ep, 1) = Location;
%     DoGMaximum.Scale(sp:ep, 1) = ScaleParameter(j) + PointCloud.ResolutionLocal(Location,1);
%     DoGMaximum.Level(sp:ep, 1) = j;
%     
%     
%     waitbar(j/NumDoGLevels, WaitBar, sprintf('Checking DoG %i of %i', j, NumDoGLevels));
% end
% 
% waitbar(j/NumDoGLevels, WaitBar, sprintf('Checking DoG for Keypoint Complete'));
% close(WaitBar)



% 
% for i = 1 : NumVertices
%     
%     
%     CurrentNeighbors = Neighbors{i,1};
%     
%     
%     
%     CurrentValue = DoG(i,j);
%     
%     SurroundingValuesCurrent = DoG(CurrentNeighbors,j);
%     
%     
%     if all(CurrentValue < SurroundingValuesCurrent)
%         isDoGMax = 1;
%     end
%     
%     
%     
%     if isDoGMax
%         NumDoGMax = NumDoGMax + 1;
%         DoGMaximum.Scale(NumDoGMax, 1) = ScaleParameter(j) + PointCloud.ResolutionLocal(i,1);
%         DoGMaximum.Location(NumDoGMax, 1) = i;
%         DoGMaximum.Level(NumDoGMax, 1) = j;
%         DoGMaximum.LocationCell{j,1} = [DoGMaximum.LocationCell{j,1}; i];
%         isDoGMax = 0;
%     end
%     
%     
% end
% 
% 
% end




% Remove empty elements from initilization

% ZeroLogic = (DoGMaximum.Location == 0);
% DoGMaximum.Scale(ZeroLogic,:) = [];
% DoGMaximum.Location(ZeroLogic,:) = [];
% DoGMaximum.Level(ZeroLogic,:) = [];


end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DoG

NumLevels = size(Signal,2);

DoG = zeros(PointCloud.Count, NumLevels-1);

for i = 1 : NumLevels - 1
    
    DoG(:,i) = Signal(:,i+1)- Signal(:,i);
%     DoG(:,i) = (Signal(:,i+1)- Signal(:,i)) * ScaleParameter(i)^2 / (ScaleParameter(i+1)^2 - ScaleParameter(i)^2);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keypoint


MaxNumberFeatures = round(NumLevels * PointCloud.LocationCount);
% Initialize space for the keypoints
Keypoint.Scale = zeros(MaxNumberFeatures,1);
Keypoint.Location = zeros(MaxNumberFeatures,1);
Keypoint.Level = zeros(MaxNumberFeatures,1);
Keypoint.Sign = zeros(MaxNumberFeatures,1);

NumFeatures = 0;
isFeaturePoint = 0;
WaitBar = waitbar(0, sprintf('Checking Vertex %i of %i', 0, PointCloud.LocationCount));

for i = 1 : PointCloud.LocationCount
    CurrentNeighbors = Neighbors.Distance{i,1};
    
    for j = 2 : NumLevels - 2
        MHEdge = 0;
        
        CurrentValue = DoG(i,j);
        
%         SurroundingValuesBelow = DoG([i,CurrentNeighbors],j-1);
        SurroundingValuesBelow = DoG(i,j-1);
        SurroundingValuesCurrent = DoG(CurrentNeighbors,j);
%         SurroundingValuesAbove = DoG([i,CurrentNeighbors],j+1);
        SurroundingValuesAbove = DoG(i,j+1);
        
        AllSign = sign([CurrentValue; SurroundingValuesCurrent; DoG(i, j+1); DoG(i, j-1)]');
        SignProd = movprod(AllSign,2,2);
        
        if nnz(SignProd(2:end) < 0) >= 2
            MHEdge = 1;
        end
        
        
        if ~MHEdge
            if CurrentValue < SurroundingValuesBelow
                if CurrentValue < SurroundingValuesAbove
                    if all(CurrentValue < SurroundingValuesCurrent)
                    isFeaturePoint = 1;
                    isLess = 1;
                    end
                end
            elseif CurrentValue > SurroundingValuesBelow
                if CurrentValue > SurroundingValuesAbove
                    if all(CurrentValue > SurroundingValuesCurrent)
                        isFeaturePoint = 1;
                        isMore = 1;
                    end
                end
            end
        end
        
        
        if isFeaturePoint
            NumFeatures = NumFeatures + 1;
            Keypoint.Scale(NumFeatures, 1) = ScaleParameter(j) + PointCloud.ResolutionLocal(i,1);
            Keypoint.Location(NumFeatures, 1) = i;
            Keypoint.Level(NumFeatures, 1) = j;
            isFeaturePoint = 0;
            if isLess
                Keypoint.Sign(NumFeatures,1) = -1;
            elseif isMore
                Keypoint.Sign(NumFeatures,1) = 1;
            end
            isLess = 0;
            isMore = 0;
        end
        
        
    end
    waitbar(i/PointCloud.LocationCount, WaitBar, sprintf('Checking Vertex %i of %i', i, PointCloud.LocationCount));
    
end

waitbar(i/PointCloud.LocationCount, WaitBar, sprintf('Checking Vertex for Keypoint Complete'));
close(WaitBar)


% Remove empty elements from initilization
ZeroLogic = (Keypoint.Location == 0);

Keypoint.Scale(ZeroLogic,:) = [];
Keypoint.Location(ZeroLogic,:) = [];
Keypoint.Level(ZeroLogic,:) = [];
Keypoint.Sign(ZeroLogic,:) = [];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



DoGLenth = size(DoG,2);

[MaxDoGVal, MaxDoGInd] = max(DoG,[],2);
[MinDoGVal, MinDoGInd] = min(DoG,[],2);

% [MaxDoGVal, MaxDoGInd] = max(Signal,[],2);
% [MinDoGVal, MinDoGInd] = min(Signal,[],2);


LocationsMax = [1:PointCloud.LocationCount]';
LocationsMin = [1:PointCloud.LocationCount]';

OutOfBoundsMaxLogic = (MaxDoGInd == 1) | (MaxDoGInd == DoGLenth);
OutOfBoundsMinLogic = (MinDoGInd == 1) | (MinDoGInd == DoGLenth);

OutOfBoundsMaxLowLogic = MaxDoGInd <= 85;
OutOfBoundsMinLowLogic = MinDoGInd <= 85;

LocationsMax(OutOfBoundsMaxLogic | OutOfBoundsMaxLowLogic) = [];
LocationsMin(OutOfBoundsMinLogic | OutOfBoundsMinLowLogic) = [];

MaxDoGVal(OutOfBoundsMaxLogic | OutOfBoundsMaxLowLogic) = [];
MaxDoGInd(OutOfBoundsMaxLogic | OutOfBoundsMaxLowLogic) = [];

MinDoGVal(OutOfBoundsMinLogic | OutOfBoundsMinLowLogic) = [];
MinDoGInd(OutOfBoundsMinLogic | OutOfBoundsMinLowLogic) = [];


QuantMax = quantile(MaxDoGVal, [0.25,0.5,0.75]);
QuantMin = quantile(MinDoGVal, [0.25,0.5,0.75]);


RemoveMin = [];
for i = 1 : length(LocationsMin)
    
    CurrentNeighbors = Neighbors.Distance{i,1};
    CurrentDoGVal = MinDoGVal(i);
    CurrentDoGInd = MinDoGInd(i);
    
    AllSign = sign([CurrentDoGVal;DoG(CurrentNeighbors, CurrentDoGInd);DoG(LocationsMin(i), CurrentDoGInd+1);DoG(LocationsMin(i), CurrentDoGInd-1)]');
    
    SignProd = movprod(AllSign,2,2);
    
    if nnz(SignProd(2:end) < 0) >= 2
        RemoveMin = [RemoveMin,i];
    end
       
   
    if ~all( CurrentDoGVal < DoG(CurrentNeighbors, CurrentDoGInd) ) %|| ~(CurrentDoGVal < DoG(LocationsMin(i), CurrentDoGInd+1)) || ~(CurrentDoGVal < DoG(LocationsMin(i), CurrentDoGInd-1))
        RemoveMin = [RemoveMin,i];
    end
    
end
RemoveMin = unique(RemoveMin);
MinDoGVal(RemoveMin) = [];
MinDoGInd(RemoveMin) = [];
size(LocationsMin)
LocationsMin(RemoveMin) = [];
size(LocationsMin)

RemoveMax = [];
for i = 1 : length(LocationsMax)
    
    CurrentNeighbors = Neighbors.Distance{i,1};
    CurrentDoGVal = MaxDoGVal(i);
    CurrentDoGInd = MaxDoGInd(i);
    
	AllSign = sign([CurrentDoGVal;DoG(CurrentNeighbors, CurrentDoGInd);DoG(LocationsMax(i), CurrentDoGInd+1);DoG(LocationsMax(i), CurrentDoGInd-1)]');
    
    SignProd = movprod(AllSign,2,2);
    
    if nnz(SignProd(2:end) < 0) >= 2
        RemoveMax = [RemoveMax,i];
    end
    
    
    if ~all( CurrentDoGVal > DoG(CurrentNeighbors, CurrentDoGInd) ) %|| ~(CurrentDoGVal > DoG(LocationsMax(i), CurrentDoGInd+1)) || ~(CurrentDoGVal > DoG(LocationsMax(i), CurrentDoGInd-1))
        RemoveMax = [RemoveMax,i];
    end
    
end
RemoveMax = unique(RemoveMax);
MaxDoGInd(RemoveMax) = [];
MaxDoGVal(RemoveMax) = [];
size(LocationsMax)
LocationsMax(RemoveMax) = [];
size(LocationsMax)


figure
subplot(1,2,1)
hist(MaxDoGInd,3000)
title('Max')
subplot(1,2,2)
hist(MinDoGInd,3000)
title('Min')


% figure
% subplot(1,2,1)
% boxplot(MaxDoGVal)
% title('Max')
% subplot(1,2,2)
% boxplot(MinDoGVal)
% title('Min')


% KpMin = MinDoGVal < (QuantMin(3) - 1.5*(QuantMin(3) - QuantMin(1)));
% KpMax = MaxDoGVal > (QuantMax(3) + 1.5*(QuantMax(3) - QuantMax(1)));
% 
% KpMin = MinDoGVal > QuantMin(2);
% KpMax = MaxDoGVal < QuantMax(2);

% KpMin = MinDoGVal < QuantMin(3);
% KpMax = MaxDoGVal > QuantMax(3);
% 
% MaxDoGVal(~KpMax) = [];
% MaxDoGInd(~KpMax) = [];
% 
% MinDoGVal(~KpMin) = [];
% MinDoGInd(~KpMin) = [];
% 
% LocationsMax(~KpMax) = [];
% LocationsMin(~KpMin) = [];


% size(LocationsMax)
% size(LocationsMin)
% 
% size(MaxDoGVal)
% size(MinDoGVal)



fig=figure('Units', 'Normalized', 'OuterPosition', [0 0 0.98 0.98]);

for i = 1 : NumLevels - 1
    
    trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), Signal(:,i), 'EdgeColor', 'none');
    % trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), DoG(:,10), 'EdgeColor', 'none');
    view(180, -90)
    ax = gca;
    axis equal
    axis off
    hold on
    cbar = colorbar('FontSize',15);
    title(sprintf('Diffusion Step: %d, Absolute Scale: %0.2f',i, ScaleParameterAbsolute(i+1)),'fontsize',25)
    
%     CurrentLocationsMin = LocationsMin(MinDoGInd==i);
%     CurrentLocationsMax = LocationsMax(MaxDoGInd==i);
    
    CurrentLocationsMin = Keypoint.Location((Keypoint.Level==i) & (Keypoint.Sign<0));
    CurrentLocationsMax = Keypoint.Location((Keypoint.Level==i) & (Keypoint.Sign>0));
    
    plot3(PointCloud.Location(CurrentLocationsMin,1), PointCloud.Location(CurrentLocationsMin,2),...
        PointCloud.Location(CurrentLocationsMin,3), 'b.','MarkerSize', 10);
    
    
    
    plot3(PointCloud.Location(CurrentLocationsMax,1), PointCloud.Location(CurrentLocationsMax,2),...
        PointCloud.Location(CurrentLocationsMax,3), 'r.','MarkerSize', 10);
    
    
    ylabel(cbar, 'Mean Curvature', 'FontSize',30);
    caxis([-0.15,0.2])
    
    drawnow
    
    clf(fig);
    
end


figure
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), Signal(:,1), 'EdgeColor', 'none');
% trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), DoG(:,10), 'EdgeColor', 'none');
view(180, -90)
ax = gca;
axis equal
axis off
hold on
cbar = colorbar('FontSize',15);
title(sprintf('Diffusion Step: %d, Absolute Scale: %0.2f',i, ScaleParameterAbsolute(i+1)),'fontsize',25)


plot3(PointCloud.Location(LocationsMin,1), PointCloud.Location(LocationsMin,2),...
    PointCloud.Location(LocationsMin,3), 'b.','MarkerSize', 10);



plot3(PointCloud.Location(LocationsMax,1), PointCloud.Location(LocationsMax,2),...
    PointCloud.Location(LocationsMax,3), 'r.','MarkerSize', 10);

















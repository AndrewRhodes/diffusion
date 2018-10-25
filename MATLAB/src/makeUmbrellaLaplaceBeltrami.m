






function [ItL, LBM] = makeUmbrellaLaplaceBeltrami(PointCloud, Neighbors, StepSize, Alpha, OrderBDF)


% Diagonal weights and locations
DiagWeights = - ones(PointCloud.LocationCount,1);
DiagLocation_i = (1:PointCloud.LocationCount)';
DiagLocation_j = (1:PointCloud.LocationCount)';


% Anonymous Functions for making weights
makelocations = @(index, neigh) index * ones(length(neigh),1);
makeweights = @(neigh) ones(length(neigh),1) / length(neigh);

% Off diagonal weights and locations
OffDiagLocations_i = cellfun(makelocations, num2cell(1:PointCloud.LocationCount)', Neighbors, 'Uni', 0);
OffDiagLocations_j = cellfun(@transpose,Neighbors,'Uni',0);
OffDiagWeights = cellfun(makeweights, Neighbors, 'Uni', 0);


Loc_i = [DiagLocation_i; vertcat(OffDiagLocations_i{:})];
Loc_j = [DiagLocation_j; vertcat(OffDiagLocations_j{:})];
Weights = [DiagWeights; vertcat(OffDiagWeights{:})];



% Construct the LBM
LBM = sparse(Loc_i, Loc_j, Weights, PointCloud.LocationCount, PointCloud.LocationCount, length(Weights));



ItL = cell(OrderBDF,1);

ItL{1,1} = speye(PointCloud.LocationCount, PointCloud.LocationCount) - Alpha * StepSize * LBM;


if OrderBDF == 2 || OrderBDF == 3 || OrderBDF == 4
%     I23tL
    ItL{2,1} = speye(PointCloud.LocationCount, PointCloud.LocationCount) - (2/3) * Alpha * StepSize * LBM;
end

if OrderBDF == 3 || OrderBDF == 4
%     I611tL
    ItL{3,1} = speye(PointCloud.LocationCount, PointCloud.LocationCount) - (6/11) * Alpha * StepSize * LBM;
end

if OrderBDF == 4
%     I1225tL
    ItL{4,1} = speye(PointCloud.LocationCount, PointCloud.LocationCount) - (12/25) * Alpha * StepSize * LBM;
end





end

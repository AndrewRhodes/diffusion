



function PointCloud = findLocalResolution(PointCloud, Neighbors)


PointCloud.ResolutionLocal = zeros(PointCloud.LocationCount,1);


for i = 1 : PointCloud.LocationCount
    
    demean = bsxfun(@minus, PointCloud.Location(Neighbors{i,1},:), PointCloud.Location(i,:));
    
    PointCloud.ResolutionLocal(i,1) = sum(sqrt(sum(demean.^2,2)),1) / size(demean,1);
    
end




end
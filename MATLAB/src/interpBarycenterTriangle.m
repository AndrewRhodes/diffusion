% Andrew Rhodes
% ASEL 2017
% interpolation across a triangular face using barycentric coordinates
%
% input[PointCloud]: struct of Point Cloud. Must contain fields 'FaceArea' 
% with size[fx3], 'Location' with size [nx3], and 'Face' with size[fx3]
% 
% input[ClosestPoint]: size[mx3] clostest points on the object surface
% from the ebedded grid.
% 
% input[cpFace]: size[mx1] face upon which the closest points reside
% 
% ouput[FaceInterpolateWeights]: size[mxn] where m is the number of
% closestpoints and n is the number of vertices on the mesh




function FaceInterpolateWeights = interpBarycenterTriangle(PointCloud, ClosestPoint, cpFace)


if ~isfield(PointCloud, 'FaceArea')
    warning('PointCloud needs the field ''FaceArea''. Calling function findFaceAreas.')
    PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
end
if ~isfield(PointCloud, 'Location')
    error('PointCloud needs the field ''Location''.')
end
if ~isfield(PointCloud, 'Face')
    error('PointCloud needs the field ''Face''.')
end


% Seperate faces
F1 = PointCloud.Face(cpFace,1);
F2 = PointCloud.Face(cpFace,2);
F3 = PointCloud.Face(cpFace,3);

% Seperate vertices
V1 = PointCloud.Location(F1,:);
V2 = PointCloud.Location(F2,:);
V3 = PointCloud.Location(F3,:);


CPArea = PointCloud.FaceArea(cpFace,1);


V1CP = V1 - ClosestPoint;
V2CP = V2 - ClosestPoint;
V3CP = V3 - ClosestPoint;


% v1 triangle area
WeightV1 = sqrt( sum( cross(V2CP, V3CP).^2 ,2 ) ) ./ (2*CPArea);

% v2 triangle area
WeightV2 = sqrt( sum( cross(V1CP, V3CP).^2 ,2 ) ) ./ (2*CPArea);

% v3 triangle area
WeightV3 = sqrt( sum( cross(V1CP, V2CP).^2 ,2 ) ) ./ (2*CPArea);

cpLength = length( cpFace );

cpVector = ( 1 : cpLength )';


FaceInterpolateWeights = sparse( [cpVector; cpVector; cpVector], [F1; F2; F3], [WeightV1; WeightV2; WeightV3] , cpLength, PointCloud.LocationCount);


% Sanity Check
nnzFaceInterpolateWeights = nnz(FaceInterpolateWeights);
nnzWeight = nnz(WeightV1) + nnz(WeightV2) + nnz(WeightV3);

if nnzFaceInterpolateWeights ~= nnzWeight
    warning('Triagle Barycentric Interpolation has missing weights.')
end


end

















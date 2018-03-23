% Andrew Rhodes
% ASEL 2017
% Find the normal vector of each face
%
% input[vertices]: size[nx3]
% input[Faces]: size[mx3]
% 
% ouput[FaceNormals]: unit vectors. size[mx3]



function FaceNormals = findFaceNormals(vertices, faces)


% Seperate faces
F1 = faces(:,1);
F2 = faces(:,2);
F3 = faces(:,3);

% Seperate vertices
V1 = vertices(F1,:);
V2 = vertices(F2,:);
V3 = vertices(F3,:);


% Edge Vectors
E1 = V1-V2;
% E2 = V2-V3;
E3 = V3-V1;

% Normalize Edge Vectors
E1 = bsxfun(@rdivide, E1, sqrt(sum(E1.^2,2)));
% E2 = bsxfun(@rdivide, E2, sqrt(sum(E2.^2,2)));
E3 = bsxfun(@rdivide, E3, sqrt(sum(E3.^2,2)));

% Calculate Face Normals
FaceNormals = cross(E1, E3);

% Normalize Face Normals
FaceNormals = bsxfun(@rdivide, FaceNormals, sqrt(sum(FaceNormals.^2,2)));



end
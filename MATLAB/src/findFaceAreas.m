% Andrew Rhodes
% ASEL 2017
% Find the area of each face
%
% input[vertices]: size[nx3]
% input[Faces]: size[mx3]
% 
% ouput[FaceArea]: unit vectors. size[mx3]





function FaceArea = findFaceAreas(vertices, faces)


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


% Calculate Face Normals
FaceArea = sqrt(sum(cross(E1, E3).^2,2));





end
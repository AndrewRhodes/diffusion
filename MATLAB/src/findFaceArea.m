% Andrew Rhodes
% ASEL 2017
% Find the area of the faces
%
% input[vertices]: size[mx3]
% input[faces]: size[nx3]
%
% output[FaceArea]: The area of the faces. Same size as input[faces].
%                   size[nx1]



function FaceArea = findFaceArea(vertices, faces)


% Ensure that vertices and faces are row major 
[m,n] = size(vertices);
if m < n
    vertices = vertices';
end

[m,n] = size(faces);
if m < n
    faces = faces';
end



V1 = vertices(faces(:,1),:);
V2 = vertices(faces(:,2),:);
V3 = vertices(faces(:,3),:);

FaceArea = sqrt( sum(bsxfun(@cross, (V2-V1), (V3-V1)).^2, 2)) / 2;


end



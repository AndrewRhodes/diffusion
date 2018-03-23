% Andrew Rhodes
% ASEL 2017
% find the adjacent neighbors for each point
%
% input[faces]: size[mx1]
% input[vertices]: size[nx1]
% 
% output[Neighbors]: cell size[nx1]



function Neighbors = findAdjacentNeighbors(faces, vertices)


NumVertices = size(vertices, 1);
NumFaces = size(faces,1);

%% Option 1: For loop

Neighbors = cell(NumVertices,1);

for i = 1 : NumFaces
    Neighbors{faces(i,1)} = [Neighbors{faces(i,1)} [faces(i,2) faces(i,3)]];
    Neighbors{faces(i,2)} = [Neighbors{faces(i,2)} [faces(i,3) faces(i,1)]];
    Neighbors{faces(i,3)} = [Neighbors{faces(i,3)} [faces(i,1) faces(i,2)]];
end

Neighbors = cellfun(@unique, Neighbors, 'UniformOutput', 0);


end
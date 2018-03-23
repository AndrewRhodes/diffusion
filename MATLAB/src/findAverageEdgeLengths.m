% Andrew Rhodes
% ASEL 2017
% Find the edge lengths
%
% input[vertices]: object vertices. size[mx3]
% input[faces]: object face matrix. size[nx3]
% input[type]: Two types available. Default[Type 1] (Recommended)
%               Type 1: Find the average edge length of the entire model. 
%                       This would result in isotropic smoothing.
%                       
%               Type 2: Find the average edge length of each face.
%                       This results in anisotropic adaptive smoothing.
% 
% output[AverageEdgeLength]: The average edge length of (Type 1) all the 
%                            faces as a scalar or (Type 2) each face as a 
%                            [nx1] matrix of scalars.
%


function AverageEdgeLength = findAverageEdgeLengths(vertices, faces, type)


% Ensure that vertices and faces are row major 
[m,n] = size(vertices);
if m < n
    vertices = vertices';
end

[m,n] = size(faces);
if m < n
    faces = faces';
end

NumFaces = size(faces,1);

if type == 1 % Isotropic
    
    % Find the average of all the edge lengths
    
    % Mean Edge Length
    % Used by Belkin
%     EdgeLength = (1/(3*NumFaces)) *...
%                  sum(sqrt(sum(([vertices(faces(:,1),:) - vertices(faces(:,2),:);
%                  vertices(faces(:,2),:) - vertices(faces(:,3),:);
%                  vertices(faces(:,3),:) - vertices(faces(:,1),:)]).^2, 2)), 1);
    
    % Median Edge Length
    % I prefer this method
    % Or, disregard distances outside of the 3*IQR and then find mean.
     EdgeLength = median(sqrt(sum(([vertices(faces(:,1),:) - vertices(faces(:,2),:);
                  vertices(faces(:,2),:) - vertices(faces(:,3),:);
                  vertices(faces(:,3),:) - vertices(faces(:,1),:)]).^2, 2)));
    
    
elseif type == 2 % Anisotropic
    
    % Mean Edge Length
    % Used by Belkin
    EdgeLength = sum([sqrt(sum((vertices(faces(:,1),:) - vertices(faces(:,2),:)).^2,2)),...
                 sqrt(sum((vertices(faces(:,2),:) - vertices(faces(:,3),:)).^2,2)),...
                 sqrt(sum((vertices(faces(:,3),:) - vertices(faces(:,1),:)).^2,2))],2)/3;
             
             
    % Median Edge Length
    % I prefer this method
%     EdgeLength = median([sqrt(sum((vertices(faces(:,1),:) - vertices(faces(:,2),:)).^2,2)),...
%                  sqrt(sum((vertices(faces(:,2),:) - vertices(faces(:,3),:)).^2,2)),...
%                  sqrt(sum((vertices(faces(:,3),:) - vertices(faces(:,1),:)).^2,2))],2);
    
             
             
else
    error('Please enter a valid type ''1'' or type ''2''.')
end


AverageEdgeLength = EdgeLength;


end












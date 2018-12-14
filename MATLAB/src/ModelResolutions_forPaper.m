% Andrew Rhodes
% Dec. 2018

% Find the model resolutions for placement in table
% Table 1 in journal paper


clc

Models = {'armadillo', 'bunny', 'buddha', 'dragon', 'itokawa'};
mod = {'Armadillo_e1', 'Bunny_e1', 'Buddha_e1', 'Dragon_e1', 'Itokawa_e1'};

n(:,1) = [103782; 121080; 138378; 155674; 172972; 190268; 207566; 224864]; %armadillo
n(:,2) = [57643; 59033; 61117; 62505; 64589; 65977; 68062; 69451]; %bunny
n(:,3) = [48737; 51985; 54153; 57409; 59577; 62819; 64999; 70430]; % buddha
n(:,4) = [47740; 50353; 52089; 54690; 56429; 59025; 60759; 62497]; % dragon
n(:,5) = [94366; 98302; 102234; 108132; 114030; 117962; 121896; 127794]; %itokawa

Resolutions = zeros(size(n,1),length(Models));

for j = 1 : length(Models)
    
    for i = 1 : size(n,1)
        
        sprintf('%s, Face: %d', Models{j}, n(i,j))
        
        Name = strcat(ProjectRoot,'/models/object/',Models{j},'/Sample/',...
            mod{j},'_',num2str(n(i,j)),'.ply');
        
        [PointCloud.Location, PointCloud.Face, PointCloud.Normal, PointCloud.Signal]...
            = read_ply_all_elements( Name );
        
        PointCloud.LocationCount = size(PointCloud.Location,1);
        PointCloud.FaceCount = size(PointCloud.Face, 1);
        PointCloud = findMeshResolution(PointCloud, 'Model');
        
        Resolutions(i,j) = PointCloud.Resolution;
        
    end
end





global ProjectRoot; % Additional Paths


modelData = {'bunny','armadillo','buddha','dragon','itokawa'}; %
ModelNames = {'Bunny_e1', 'Armadillo_e1_100000', 'Buddha_e1_50000', 'Dragon_e1_50000','Itokawa_e1_80000'}; %

NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];

for i = 1 : length(modelData)
    
    Model = strcat(modelData{i},'/',ModelNames{i})
    
    FileLocationModel = strcat(ProjectRoot,'/models/object/');
    FileNameModelPly = strcat(Model,'.ply');
    FileNameModelOff = strcat(Model,'.off');
    
    [PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
    
    PointCloud.LocationCount = size(PointCloud.Location,1);
    PointCloud.FaceCount = size(PointCloud.Face, 1);
    PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
    PointCloud = findMeshResolution(PointCloud, 'Model');
    PointCloud = findMeshNormals(PointCloud);
    
    
    load(strcat(ModelNames{i},'_Neighbors.mat'),'Neighbors')
    PointCloud = findLocalResolution(PointCloud, Neighbors.Connect)
    
    
    for j = 1 : length(NoiseVec)
    
        Location = PointCloud.Location + NoiseVec(j) * PointCloud.Resolution * ( randn(PointCloud.LocationCount,1) .* PointCloud.Normal );
        File = strcat(ProjectRoot, '/models/object/',modelData{i},'/',ModelNames{i},'_sigma',num2str(j));
        
        ply_write(strcat(File,'.ply'), PointCloud.Face, Location)
        save_off(Location, PointCloud.Face, strcat(File,'.off'))
        
    end
    
    clear PointCloud Neighbors
end









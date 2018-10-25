


modelData = {'bunny','armadillo','buddha','dragon','itokawa'}; %
ModelNames = {'Bunny_e1', 'Armadillo_e1_100000', 'Buddha_e1_50000', 'Dragon_e1_50000','Itokawa_e1_80000'}; %
ModelVertices = [34834, 50002, 24939, 24956, 39998];
NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];


global ProjectRoot;


for i = 1 : length(ModelNames)
    
    for j = 1 : length(NoiseVec)
        clear Curvature
        
        Path = strcat(ProjectRoot,'/models/object/',modelData{i},'/',ModelNames{i},'_sigma',num2str(j),'_MK.ply');
        
        fid = fopen(Path,'r+')
        
        c = 1;
        
        while c
            if strcmpi(fgetl(fid),'end_header')
                c = 0;
            end
        end
        
        
        Curvature = zeros(ModelVertices(i),1);
        
        for k = 1 : ModelVertices(i)
            
            buf = fgetl(fid);
            buf = strsplit(buf);
            Curvature(k,1) = str2double(buf{4});
            
        end
        
        fclose(fid);
        
        save(strcat(ProjectRoot,'/main/DE/keypointdata/',modelData{i},'/',ModelNames{i},'_Curvature_sigma',num2str(j),'.mat'), 'Curvature', '-v7.3')
        
        
    end
end















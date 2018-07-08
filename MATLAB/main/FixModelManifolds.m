




PointCloudNew = PointCloud;
RemoveIndex=[];
for i = 1 : PointCloud.FaceCount
    
    CurrentFace = PointCloud.Face(i,:);
    
    MatchingFaces = (PointCloud.Face == CurrentFace(1)) ...
        | (PointCloud.Face == CurrentFace(2)) ...
        | (PointCloud.Face == CurrentFace(3));
    
    LogicMatchingFaces = sum(MatchingFaces,2) == 3;
    nnzLogicMatchingFaces = nnz(LogicMatchingFaces);
    
    if nnzLogicMatchingFaces > 1
        i
        nnzLogicMatchingFaces
        IndexMatchingFace = find(LogicMatchingFaces);
        RemoveIndex = [RemoveIndex; IndexMatchingFace(2:end)]
    end
    
end


RemoveIndex
RemoveIndex = unique(RemoveIndex)


PointCloudNew.Face(RemoveIndex,:) = [];


ModelName = 'Buddha_e1';
DownSampleFacesNum = 100000;



ply_write(strcat(ModelName,'_',num2str(DownSampleFacesNum),'.ply'), PointCloud.Face, PointCloud.Location)
save_off(PointCloud.Location, PointCloud.Face, strcat(ModelName,'_',num2str(DownSampleFacesNum),'.off'))


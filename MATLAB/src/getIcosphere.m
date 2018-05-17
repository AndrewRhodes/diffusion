

function PointCloud = getIcosphere(FileName, NumberDivisions)



if exist(FileName, 'file')
    
    [PointCloud.Location, PointCloud.Face] = read_off( FileName );

    [m, n] = size(PointCloud.Location);
    if m < n
        PointCloud.Location = PointCloud.Location';
    end
    
    [m, n] = size(PointCloud.Face);
    if m < n
        PointCloud.Face = PointCloud.Face';
    end
    

else
    
    [Location, Faces] = icosphere(NumberDivisions);
    [VerticesOut, FacesOut] = clearMeshDuplicates(Location, Faces );
    
    save_off(VerticesOut, FacesOut, fullfile( FileName ) )
    
    PointCloud.Face = FacesOut;
    PointCloud.Location = VerticesOut;
    
    
end



end


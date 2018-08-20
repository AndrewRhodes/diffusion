



global ProjectRoot;


Path = strcat(ProjectRoot,'/models/object/itokawa/Itokawa_e1_100000_MK.ply');

fid = fopen(Path,'r+')

c = 1;

while c
    if strcmpi(fgetl(fid),'end_header')
        c = 0;
    end
end

NumVert = 50001;

Curvature = zeros(NumVert,1);

for i = 1 : NumVert
    
    buf = fgetl(fid);
    buf = strsplit(buf);
    Curvature(i,1) = str2double(buf{4});

end

save ItokawaCurvature_e1_100000 Curvature


fclose(fid);

















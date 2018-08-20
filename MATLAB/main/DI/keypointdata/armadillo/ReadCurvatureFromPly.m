



fid = fopen(fullfile( FileLocationModel, FileNameModelPly ));
c = 1;

while c
    Buf = fscanf(fid,'%s',1);
%     pause
    if strcmpi(Buf, 'end_header')
        c=0;
        r=1;
        fgetl(fid);
    end
end
Curvature = zeros(50002,1);

for i = 1 : 50002
    CurrentLine = fgetl(fid);
    CurrentLineSep = strsplit(CurrentLine);
    Curvature(i,1) = str2double(CurrentLineSep{1,4});
%     pause
end

fclose(fid);
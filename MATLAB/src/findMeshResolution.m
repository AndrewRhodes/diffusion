% Andrew Rhodes
% ASEL
% November 2015
% -----
% input[PointCloudIn]:
% input[Type]: string 'Scene' or 'Model'
% output[PointCloudIn]: With new field .Resolution

function PointCloudIn = findMeshResolution(PointCloudIn, Type)

if strcmp(Type, 'Model')
    Fields = fieldnames(PointCloudIn);
    
    if ~strcmp(Fields, 'Location')
        error('PointCloudIn must have a data field ''Location''')
    elseif ~strcmp(Fields, 'Face')
        error('PointCloudIn must have a data field ''Face''')
    end

    A(:,:,1) = PointCloudIn.Location(PointCloudIn.Face(:,1),:);
    
    for i = 1 : 3
        j = mod(i,3) + 1;
        A(:,:,i+1) = PointCloudIn.Location(PointCloudIn.Face(:,j),:);
    end
    
    EdgeLengths = reshape(squeeze(sqrt(sum(diff(A,1,3).^2,2))),[],1);
    
    % As defined by A.J. Spin Images
    PointCloudIn.Resolution = median(EdgeLengths);
elseif strcmp(Type, 'Scene')
    
    % Use interquantile range to detect outliers. Thanks for the idea A.J.
    EdgeLengths = sqrt(sum(diff(PointCloudIn.Location,1).^2,2));
    Qs = quantile(EdgeLengths, [0.25 0.5 0.75]); % Find quantiles
    % Find the median of edge lenghts that are not above the 
    % 75% quantile + 3*(interquantile range)
    PointCloudIn.Resolution = median(EdgeLengths(EdgeLengths<= (Qs(3)+3*(Qs(3)-Qs(1))) ));
end
end
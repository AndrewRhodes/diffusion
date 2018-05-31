

function [ClosestPoint, Distance] = findClosestPoint(PointCloudIn, Point, FaceIndices)

NumFaces = length(FaceIndices);
dist = zeros(1,NumFaces);
cp = zeros(NumFaces,3);

for i = 1 : NumFaces
    
    % Collect vertices of Face
    v1 = PointCloudIn.Location(PointCloudIn.Face(FaceIndices(i)),:);
    v2 = PointCloudIn.Location(PointCloudIn.Face(FaceIndices(i)),:);
    v3 = PointCloudIn.Location(PointCloudIn.Face(FaceIndices(i)),:);
    
    % Move v1 to the origin
    
    v1p = Point - v1; % a
    v2p = Point - v2; % q
    v3p = Point - v3; % r
    
    % Evaluate matrix entries
    a12 = v2p * v3p';
    A = [(v2p * v2p'), a12; a12, (v3p * v3p')];
    b= [v1p * v2p'; v1p * v3p'];
    
    x = A \ b;
    
    cp(i,:) = [(x(1) * v2p(1)) + (x(2) * v3p(1)), ...
        (x(1) * v2p(2)) + (x(2) * v3p(2)), ...
        (x(1) * v2p(3)) + (x(2) * v3p(3))];
    
    
    if (x(1)<0) && (x(2)<0) && (sum(x)<=1)
        cp(i,:) = zeros(1,3);
    elseif (x(1)>=0) && (x(2)<0) && (sum(x)<=1)
        cp(i,:) = ProjectOnSegment(cp(i,:), zeros(1,3), v2p);
    elseif (x(1)>=0) && (x(2)<0) && (sum(x)>1)
        cp(i,:) = v3p;
    elseif (x(1)<0) && (x(2)>=0) && (sum(x)<=1)
        cp(i,:) = ProjectOnSegment(cp(i,:),v3p,zeros(1,3));
    elseif (x(1)>=0) && (x(2)>=0) && (sum(x)<=1)
%         warning("Shoudn't be here");
    else
%         warning("Didn't account for something");
    end
    
    if abs(1/det(A)) > 1e15
        dist(i) = 10000;
        warning("Panic!");
    else
        dist(i) = sqrt(sum((v1p-cp(i,:)).^2));
    end
    
    % Shift everything back
    cp(i,:) = cp(i,:) + v1;

end


[Distance, Index] = min(dist);


ClosestPoint = cp(Index,:);



end

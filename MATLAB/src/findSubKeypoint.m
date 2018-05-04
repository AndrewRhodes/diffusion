









function SubKeypoint = findSubKeypoint(PointCloudIn, NormalRotationsIn, DoG, Keypoint, Neighbors, )

if ~isfield(PointCloudIn, 'Location')
    error('PointCloudIn must contain the field ''Location''.')
end

if ~isfield(PointCloudIn, 'Normal')
    error('PointCloudIn must contain the field ''Normal''.')
end

if ~isfield(PointCloudIn, 'NormalRotationMatix')
    warning('Adding field to ''NormalRoationMatrix'' to PointCloudIn')
    PointCloudIn = findNormalsRotation(PointCloudIn);
end

NumKeypoints = length(Keypoint.Location);

SubKeypoint.Location = zeros(NumKeypoints,3);
SubKeypoint.Scale = zeros(NumKeypoints,1);
% Subkeypoint.Level = zeros(NumKeypoints,1);


for i = 1 : NumKeypoints
    
    CurrentVertex = PointCloudIn.Location(Keypoint.Location(i),:);
    CurrentNormal = PointCloudIn.Normal(Keypoint.Location(i),:);
    CurrentScale = Keypoint.Scale(i);
	CurrentNeighbors = Neighbors{i,1};
    
    m = 3 * length(CurrentNeighbors) - 1;
    x = zeros(m, 4);
    phi = zeros(m, 10);
    D = zeros(m, 1);
    
   for j = 1 : 3 * length(CurrentNeighbors) - 1
        x(j,1:4) = [CurrentVertex(1), CurrentVertex(2), CurrentVertex(3), CurrentScale];
        phi(j,1:10) = buildPhi(CurrentVertex(1), CurrentVertex(2), CurrentVertex(3), CurrentScale);
        D(j,1) = DoG(i, Keypoint.Level(i) );
    end
    
    PHI = [x, phi, ones(m,1), -D];
    
    [~,~,V] = svd(PHI);
    
    y = V(:,end) ./ V(end,end);
    
    [J, H] = buildJHc(y);
    
    if det(H) ~= 0
        xhattilde = -inv(H) * J;
        
        xhat = [CurrentNormal'*CurrentNormal, 0;
            0,                            1] * xhattilde;
        
        SubKeypoint.Location = xhat(1:3)';
        SubKeypoint.Scale = xhat(4);
        
    else
        
        disp(cond(H))
        error('Det(H) cannot be found.')
        
    end

    

end







end




function Phi = buildPhi(vx, vy, vz, sigma)

Phi = [vx^2, 2*vx*vy, 2*vx*vz, 2*vx*sigma, ...
       vy^2, 2*vy*vz, 2*vy*sigma, ...
       vz^2, 2*vz*sigma, ...
       sigma^2];


end


function [J, H] = buildJHc(y)

    J = y(1:4,1);
    
    H = [y(5), y(6),  y(7),  y(8);
         y(6), y(9),  y(10), y(11);
         y(7), y(10), y(12), y(13);
         y(8), y(11), y(13), y(14)];


end

















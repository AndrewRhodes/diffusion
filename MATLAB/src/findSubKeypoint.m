









function SubKeypoint = findSubKeypoint(Keypoint, ScaleParameter, DoG, PointCloud, Neighbors, NeighborFaces)

IncludeNumLevels = 5; % Must be odd
Range = floor(IncludeNumLevels/2);

if ~isfield(PointCloud, 'Location')
    error('PointCloudIn must contain the field ''Location''.')
end

% if ~isfield(PointCloudIn, 'Normal')
%     error('PointCloudIn must contain the field ''Normal''.')
% end

% if ~isfield(PointCloudIn, 'NormalRotationMatix')
%     warning('Adding field to ''NormalRoationMatrix'' to PointCloudIn')
%     PointCloudIn = findNormalsRotation(PointCloudIn);
% end

NumKeypoints = size(Keypoint.Location,1);

SubKeypoint.Location = zeros(NumKeypoints,3);
SubKeypoint.Scale = zeros(NumKeypoints,1);
SubKeypoint.Index = zeros(NumKeypoints,1);
SubKeypoint.Distance = zeros(NumKeypoints,1);
SubKeypoint.Covariance = zeros(4,4,NumKeypoints);

for i = 1 : NumKeypoints
    
    if Keypoint.Level(i) > 3
        CurrentVertex = PointCloud.Location(Keypoint.Location(i),:);
        %     CurrentNormal = PointCloudIn.Normal(Keypoint.Location(i),:);
        CurrentScale = Keypoint.Scale(i);
        CurrentNeighbors = Neighbors{i,1};
        CurrentDoG = DoG(Keypoint.Location(i), Keypoint.Level(i));
        
        m = IncludeNumLevels * length(CurrentNeighbors) + IncludeNumLevels;
        x = zeros(m, 4);
        phi = zeros(m, 10);
        D = zeros(m, 1);
        
        c = 0;
        for j = -Range : Range
            c = c +1;
            
            scale = ScaleParameter(Keypoint.Level(i)+j);
            
            vertex = CurrentVertex;
            
            x(c,1:4) = [vertex(1), vertex(2), vertex(3), scale];
            
            phi(c,1:10) = buildPhi(vertex(1), vertex(2), vertex(3), scale);
            
            D(c,1) = DoG(Keypoint.Location(i), Keypoint.Level(i)+j);
            
        end
        
        
        % Over the 3 levels: below, current, above
        for j = -Range : Range
            
            scale = ScaleParameter(Keypoint.Level(i)+j);
            
            % Over neighbors of current vertex (not including current vertex)
            for k = 1 : length(CurrentNeighbors)
                
                c = c + 1;
                
                vertex = PointCloud.Location(CurrentNeighbors(k),:);
                
                x(c,1:4) = [vertex, scale];
                
                phi(c,1:10) = buildPhi(x(c,1), x(c,2), x(c,3), x(c,4));
                
                D(c,1) = DoG(CurrentNeighbors(k), Keypoint.Level(i)+j);
                
            end
            
        end
        
        
        PHI = [x, phi, ones(m,1), -D];
        
        [~,~,V] = svd(PHI);
        
        y = V(:,end) ./ V(15,end);
        
        [J, H, c] = buildJHc(y);
        
        if abs(det(H)) > 1e-10
            
            xhattilde = - H \ J;
            
            %             xhattilde(1:3,1) = xhattilde(1:3,1) + CurrentVertex';
            %             xhattilde(4) = xhattilde(4) + CurrentScale;
            
            %%%
            % Apply cp(x) here
            if length(NeighborFaces{Keypoint.Location(i)}) >= 1
                [ClosestPoint, Distance] = findClosestPoint(PointCloud, xhattilde(1:3)', NeighborFaces{Keypoint.Location(i)});
                
                SubKeypoint.Location(i,:) = ClosestPoint;
                SubKeypoint.Scale(i,1) = xhattilde(4);
                SubKeypoint.Distance(i,1) = Distance;
                SubKeypoint.Index(i,1) = i;
                SubKeypoint.Covariance(:,:,i) = inv(H);
                SubKeypoint.DoG(i,1) = findOptimalDoG(J, H, c, [SubKeypoint.Location(i,:), SubKeypoint.Scale(i,1)]');
            end
        else
            
            warning('Det(H) cannot be found for keypoint %d', i)
            disp(cond(H))
            
        end
    end
    
    
end


end




% function [J, H, c] = buildJHc(y)
%
% J = y(1:4,1);
%
% H = [y(5), y(6),  y(7),  y(8);
%     y(6), y(9),  y(10), y(11);
%     y(7), y(10), y(12), y(13);
%     y(8), y(11), y(13), y(14)];
%
% c = y(15);
%
% end



% function Phi = buildPhi(vx, vy, vz, sigma)
%
% Phi = [vx^2, 2*vx*vy, 2*vx*vz, 2*vx*sigma, ...
%     vy^2, 2*vy*vz, 2*vy*sigma, ...
%     vz^2, 2*vz*sigma, ...
%     sigma^2];
%
%
% end



% function D_hat = findOptimalDoG(J, H, c, xhat)
%
% D_hat = c + J' * xhat + xhat' * (H \ xhat);
%
%
% end







% function [ClosestPoint, Distance] = findClosestPoint(PointCloudIn, Point, FaceIndices)
%
% NumFaces = length(FaceIndices);
% dist = zeros(1,NumFaces);
% cp = zeros(NumFaces,3);
%
% for i = 1 : NumFaces
%
%     % Collect vertices of Face
%     v1 = PointCloudIn.Location(PointCloudIn.Face(FaceIndices(i),1),:);
%     v2 = PointCloudIn.Location(PointCloudIn.Face(FaceIndices(i),2),:);
%     v3 = PointCloudIn.Location(PointCloudIn.Face(FaceIndices(i),3),:);
%
%     % Move v1 to the origin
%
%     v1p = Point - v1; % a
%     v2p = Point - v2; % q
%     v3p = Point - v3; % r
%
%     % Evaluate matrix entries
%     a12 = v2p * v3p';
%     A = [(v2p * v2p'), a12; a12, (v3p * v3p')];
%     b= [v1p * v2p'; v1p * v3p'];
%
%     x = A \ b;
%
%     cp(i,:) = [(x(1) * v2p(1)) + (x(2) * v3p(1)), ...
%         (x(1) * v2p(2)) + (x(2) * v3p(2)), ...
%         (x(1) * v2p(3)) + (x(2) * v3p(3))];
%
%
%     if (x(1)<0) && (x(2)<0) && (sum(x)<=1)
%         cp(i,:) = zeros(1,3);
%     elseif (x(1)>=0) && (x(2)<0) && (sum(x)<=1)
%         cp(i,:) = ProjectOnSegment(cp(i,:), zeros(1,3), v2p);
%     elseif (x(1)>=0) && (x(2)<0) && (sum(x)>1)
%         cp(i,:) = v3p;
%     elseif (x(1)<0) && (x(2)>=0) && (sum(x)<=1)
%         cp(i,:) = ProjectOnSegment(cp(i,:),v3p,zeros(1,3));
%     elseif (x(1)>=0) && (x(2)>=0) && (sum(x)<=1)
%         warning("Shoudn't be here");
%     else
%         warning("Didn't account for something");
%     end
%
%     if abs(1/det(A)) > 1e15
%         dist(i) = 10000;
%         warning("Panic!");
%     else
%         dist(i) = sqrt(sum((v1p-cp(i,:)).^2));
%     end
%
%
%     % Shift everything back
%     cp(i,:) = cp(i,:) + v1;
%
% end
%
%
% [Distance, Index] = min(dist);
%
%
% ClosestPoint = cp(Index,:);
%
%
%
% end





% function cpq = ProjectOnSegment(c1, p1, q1)
%
%
% cmp = c1 - p1;
% qmp = q1 - p1;
%
% lambda = (cmp * qmp') / (sum(qmp.^2));
% lambda = max(0, min(lambda, 1));
%
% cpq = p1 + lambda * qmp;
%
%
%
%
% end
%
%
%




























% Andrew Rhodes
% ASEL 2017
% Find curvature values for each vertex
% This is a re-work of 
% 'patchcurvature' written by D.Kroon University of Twente 
%
% input[PointCloudIn]: must contain fields 'Location', 'Normal'
%
% output[PK1]: Principal curvature in direction 1. size[nx1]
% output[PK2]: Principal curvature in direction 2. size[nx1]
% output[PD1]: Principal direction 1. size[nx3]
% output[PD2]: Principal direction 2. size[nx3]
% output[MK]: Mean Curvauture = 0.5*( PK1 + PK2 ). size[nx1]
% output[GK]: Gaussian Curvature = PK1*PK2. size[nx1]



function [PK1, PK2, PD1, PD2, MK, GK] = findPointCurvatures(PointCloudIn, NormalRotationsIn, Neighbors)

if ~isfield(PointCloudIn, 'Location')
    error('PointCloudIn must contain field ''Location''.')
end

if ~isfield(PointCloudIn, 'Normal')
    error('PointCloudIn must contain field ''Normal''.')
end

if ~isfield(PointCloudIn, 'LocationCount')
    PointCloudIn.LocationCount = length(PointCloud.Location);
end


PK1 = zeros(PointCloudIn.LocationCount, 1);
PK2 = zeros(PointCloudIn.LocationCount, 1);
PD1 = zeros(PointCloudIn.LocationCount, 3);
PD2 = zeros(PointCloudIn.LocationCount, 3);


for i = 1 : PointCloudIn.LocationCount
    
    % Second ring neighbors
    AllNeighbors = unique([Neighbors{Neighbors{i,1},1}]);
    
    % Third ring neighbors
%     AllNeighbors = unique([Neighbors{[Neighbors{Neighbors{i,1},1}],1}]);
    
    CurrentVertices = PointCloudIn.Location(AllNeighbors,:);
    
    RotatedVertices = CurrentVertices * NormalRotationsIn(:,:,i);
    
    % Fit quadratic patch
    % f(x,y) = ax^2 + by^2 + cxy + dx + ey + g
    f = RotatedVertices(:,1);
    x = RotatedVertices(:,2);
    y = RotatedVertices(:,3);
    
    % LLS because these are points on a model (exact)
    A = [x.^2, y.^2, x.*y, x, y, -f, ones(size(x))];
    [~,~,V] = svd(A);
    Coeff = V(:,7) / V(7,7);
    
%     A = [x.^2, y.^2, x.*y, x, y, ones(size(x))];
%     Coeff = (A'*A) \ (A'*f); % A \ f; %

    % Create Hessian
    % Second Derivative
    H = [2*Coeff(1), Coeff(3);
         Coeff(3), 2*Coeff(2)];
     
    % Find the principal curvatures and directions
    [P1, P2, K1, K2] = findPrincipalCurvature(H);
%     [vec,val] = eig(H,'vector');
    
    RotatedVec = [0, P1; 0, P2] * NormalRotationsIn(:,:,i)';
%     RotatedVec = bsxfun(@rdivide, RotatedVec, sqrt(sum(RotatedVec.^2,2)));
    
    PD1(i,:) = RotatedVec(1,:) / sqrt(sum(RotatedVec(1,:).^2));
    PD2(i,:) = RotatedVec(2,:) / sqrt(sum(RotatedVec(2,:).^2));
    PK1(i,1) = K1;
    PK2(i,1) = K2;
    

end


% Mean Curvature
MK = 0.5*(PK1+PK2);

% Gaussian Curvature
GK = PK1.*PK2;
 
end




function [D1, D2, K1, K2] = findPrincipalCurvature(H)
% Replacement of the eig function for this specific scenerio.
%
% input[H]: determinant H from which to determine eigenvalues/vectors
% output[D1]: first principal direction
% output[D2]: second principal direction
% output[K1]: first principal curvature
% output[K2]: second principal curvature

det = sqrt( (H(1)-H(4))^2 + 4*H(2)*H(3) );

lambda1 = 0.5 * ((H(1)+H(4)) + det);
lambda2 = 0.5 * ((H(1)+H(4)) - det);

vec2 = [(H(2) + H(3)), (H(4) - H(1) + det)];

magnitude = norm(vec2);

if magnitude < 1e-10
    % Rarely Occurs
    vec2 = [1, 0];
else
    vec2 = vec2 / magnitude; 
end

vec1 = [vec2(2), -vec2(1)];

% We desire that |PC1| > |PC2|
if abs(lambda1) > abs(lambda2)
    K1 = lambda1;
    K2 = lambda2;
    D1 = vec1;
    D2 = vec2;
else % abs(lambda2) > abs(lambda1)
    K1 = lambda2;
    K2 = lambda1;
    D1 = vec2;
    D2 = vec1;
end

end











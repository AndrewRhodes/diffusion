



function NormalRotations = findNormalsRotation(NormalsIn)


NumNormals = size(NormalsIn,1);




%% Find Rotation Matrix
% Rotate the normal to [1;0;0]

% Find the rotation angles and axes
RotationAxes = bsxfun(@cross, NormalsIn', [1;0,;0])';
RotationAngles = asin( min(1, sqrt( sum( RotationAxes.^2, 2) ) ) );

% Angle is too small to make a rotation
NoRotation = RotationAngles < 1e-10;
NumNoRotation = sum(NoRotation);

% Clean up the rotation axes
RotationAxes = bsxfun(@rdivide, RotationAxes, sqrt(sum(RotationAxes.^2,2)));
RotationAxes(NoRotation, :) = ones(NumNoRotation,1)*[0,0,0];


%% Option 1: Bsxfun
% Supter Fast

% Use the axes-angle to rotation matrix formula. 
% R = cos(A)*eye(3) + (1-cos(A))*V*V' - sin(A)*[Vx]
% But do all the processes in a row vector and then reshape after. 

% Skew symmetric matrix. [Vx]
SkewMat = [zeros(NumNormals,1), RotationAxes(:,3), -RotationAxes(:,2), -RotationAxes(:,3), zeros(NumNormals,1), RotationAxes(:,1), RotationAxes(:,2), -RotationAxes(:,1), zeros(NumNormals,1)];

% Sin(A)*[Vx]
SinSkew = bsxfun(@times, -sin(RotationAngles), SkewMat);

% cos(A)*eye(3)
CosEye3 = bsxfun(@times, cos(RotationAngles), reshape(eye(3),1,9));

% V*V'
OuterProd = repmat(RotationAxes,1,3) .* reshape(repmat(RotationAxes,3,1), [], 9);

% cos(A)*V*V'
CosOuterProd = bsxfun(@times, 1-cos(RotationAngles), OuterProd);

% Combine all the elements
RM = CosEye3 + CosOuterProd + SinSkew;

% Replace the non rotating points with identity
RM(NoRotation,:) = ones(NumNoRotation,1)* reshape(eye(3),1,9);

% reshape into a [3x3xn] array
NormalRotations = reshape(RM', 3, 3, NumNormals);


end

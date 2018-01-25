% Andrew Rhodes
% WVU
% Jan. 2018
%
% findDoG differences diffused signal of neighbor levels to create a
% Difference of Gaussian of the signal
% 
% input[DiffusedSignal]: the vectorized signal diffused 
% input[ScaleParameter]: list of scale parameters for each diffusion step
%
% output[DoG]: Difference of Gaussian of signals 

function DoG = findDoG(DiffusedSignal, ScaleParameter)

if length(size(DiffusedSignal)) > 2
    DiffusedSignal = reshape(DiffusedSignal, numel(DiffusedSignal(:,:,1)), size(DiffusedSignal,3));
end


if any(ismember(length(ScaleParameter), size(DiffusedSignal)))
    NumLevels = length(ScaleParameter);
else    
    error('ScaleParameter and DoG must have the same number of elements')
end


% Find the number of points
[Num1, Num2] = size(DiffusedSignal);
a = [Num1, Num2];
NumPoints = a(a~=NumLevels);



DoG = zeros(NumPoints, NumLevels-1);


for i = 1 : NumLevels-1
    
    DoG(:,i) = 2 * ( DiffusedSignal(:,i+1) - DiffusedSignal(:,i) ) * (ScaleParameter(i,1)^2 / (ScaleParameter(i+1,1)^2 - ScaleParameter(i,1)^2) );
    
end





end
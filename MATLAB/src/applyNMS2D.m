




function NMSKeypoint = applyNMS2D(Keypoint, DoG, Image, t_scale, t_range, t_rangeMethod, DoGNormalize, CompareMethod)



ImageSize = size(Image);



r = (10+1)^2 / 10;


% Find Image Gradients
[Gx,Gy] = imgradientxy(Image, 'central');
[Gxx, Gyx] = imgradientxy(Gx, 'central');
[Gxy, Gyy] = imgradientxy(Gy, 'central');


for i = 1 : Keypoint.Count
    
    % Discard points below 
    
    
    H = [Gxx(), Gxy(); Gxy(), Gyy()];
    trH = trace(H);
    detH = det(H);
    
%     needs to be true
    trH^2 / detH < r;

end





% Comet67P Only
ZeroLogic = (Image(Keypoint2D.LocationIndex) < 0.04) | ...
             (Keypoint2D.Location(:,2) <= 3) | ...
             (Keypoint2D.Location(:,1) <= 3) | ...
             (Keypoint2D.Location(:,1) >= ImageSize(1)-3) | ...
             (Keypoint2D.Location(:,2) >= ImageSize(2)-3);
         
FNames = fieldnames(Keypoint2D);
for jj = 1 : length(FNames)
    Keypoint2D.(FNames{jj})(ZeroLogic,:) = [];
end
Keypoint2D.Count = length(Keypoint2D.LocationIndex);






end
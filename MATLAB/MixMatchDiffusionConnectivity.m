



% load('SignalImplicitComet67P.mat')
load('Comet67P3DDiffusionSigma075.mat')
% load('Comet67PNeighbors3D.mat')
load('Comet67PNeighbors3DFlipped.mat')

% LaplacianScale = zeros(NumPixels, MaxLevel-1);
LaplacianScaleNormalized = zeros(NumPixels, MaxLevel-1);
% LaplaceScaleInvariant = zeros(NumPixels, MaxLevel-1);

% ScaleParamterImage = ScaleParameterSpatialImplicit;

for i = 1 : MaxLevel-1
    
    
    % Estimated Laplacian
    
    % Scale Normalized Laplacian
%     LaplacianScaleNormalized(:,i) = 2 * ( reshape(Signal2D(:,:,i+1),[],1) - reshape(Signal2D(:,:,i),[],1) ) * (ScaleParamterImage(i,1)^2 / (ScaleParamterImage(i+1,1)^2 - ScaleParamterImage(i,1)^2) );

    
     % Estimated Laplacian
    
    % Scale Normalized Laplacian
    LaplacianScaleNormalized(:,i) = 2 * ( SignalImplicit(:,i+1) - SignalImplicit(:,i) ) * (ScaleParameterSpatialImplicit(i,1)^2 / (ScaleParameterSpatialImplicit(i+1,1)^2 - ScaleParameterSpatialImplicit(i,1)^2) );

    
    
end



LaplaceDOG = LaplacianScaleNormalized;

FeaturePoint.Scale = zeros(365000,1);
FeaturePoint.Location = zeros(365000,1);

NumFeatures = 0;
for i = 1 : NumPixels
    % for i = 1 : PointCloud.LocationCount
    
    CurrentNeighbors = Neighbors3D{i,1}';
    
    %     CurrentNeighbors = Neighbors2D{i,1};
    %     CurrentNeighbors = sub2ind([MaxImageSize(1),MaxImageSize(2)], CurrentNeighbors(:,1), CurrentNeighbors(:,2));
    
    
    if Image(i) >= 0.04
        
        for j = 2 : MaxLevel - 2
            
            CurrentValue = LaplaceDOG(i,j);
            
            SurroundingValues_under = LaplaceDOG([i;CurrentNeighbors],j-1);
            SurroundingValues_same = LaplaceDOG(CurrentNeighbors,j);
            SurroundingValues_above = LaplaceDOG([i;CurrentNeighbors],j+1);
            
            
            if all(CurrentValue > SurroundingValues_same)
                
                if all(CurrentValue > SurroundingValues_under)
                    if all(CurrentValue > SurroundingValues_above)
                        NumFeatures = NumFeatures + 1;
                        % Maximum
                        FeaturePoint.Scale(NumFeatures, 1) = ScaleParameterSpatialImplicit(j,1);
                        FeaturePoint.Location(NumFeatures, 1) = i;
                    end
                end
                
            elseif all(CurrentValue < SurroundingValues_same)
                
                if all(CurrentValue < SurroundingValues_under)
                    if all(CurrentValue < SurroundingValues_above)
                        % Minimum
                        NumFeatures = NumFeatures + 1;
                        FeaturePoint.Scale(NumFeatures, 1) = ScaleParameterSpatialImplicit(j,1);
                        FeaturePoint.Location(NumFeatures, 1) = i;
                    end
                end
            end
            
            
        end
        
    end
end

FeaturePoint.Scale(FeaturePoint.Location == 0,:) = [];
FeaturePoint.Location(FeaturePoint.Location == 0) = [];


ZeroLogic = Image(FeaturePoint.Location) < 0.04;
FeaturePoint.Location(ZeroLogic,:) = [];
FeaturePoint.Scale(ZeroLogic,:) = [];

save FeaturePointComet67P3DDiffusion3DConnectivityFlipped FeaturePoint


% % 2D Diffusion, 2D Connectivity
% save FeaturePointComet67P2DDiffusion2DConnectivity FeaturePoint
% load('FeaturePointComet67P2DDiffusion2DConnectivitySigma.mat')
% FeaturePoint2Dd2Dc = FeaturePoint2D;


% 2D Diffusion, 3D Connectivity
% save FeaturePointComet67P2DDiffusion3DConnectivitySigma FeaturePoint
% load('FeaturePointComet67P2DDiffusion3DConnectivity.mat')
% FeaturePoint2Dd3Dc = FeaturePoint;


% 3D Diffusion, 2D Connectivity
% save FeaturePointComet67P3DDiffusion2DConnectivity FeaturePoint
% load('FeaturePointComet67P3DDiffusion2DConnectivity.mat')
% FeaturePoint3Dd2Dc = FeaturePoint;


% 3D Diffusion, 3D Connectivity
% save FeaturePointComet67P3DDiffusion3DConnectivitySigma FeaturePoint
% load('FeaturePointComet67P3DDiffusion3DConnectivitySigma.mat')
% FeaturePoint3Dd3Dc = FeaturePoint;

CircleAngle = linspace(0, 2*pi, 360);
CircleX = cos(CircleAngle);
CircleY = sin(CircleAngle);

FeaturePoint2Dd3Dc=FeaturePoint2Dd3Dc;
FeaturePoint3Dd3Dc=FeaturePoint3Dd3Dc;
FeaturePoint2Dd2Dc=FeaturePoint2Dd2Dc;
FeaturePoint3Dd2Dc=FeaturePoint3Dd2Dc;



clear SortValue SortOrder
[SortValue, SortOrder] = sort(FeaturePoint.Scale,'descend');


PlotFeatures = 1:100;

figure
imshow(Image)
hold on


t=0;
for i = 1 : length(PlotFeatures)
    
    [a,b] = ind2sub([MaxImageSize(1),MaxImageSize(2)], FeaturePoint.Location(SortOrder(i)));
    
    CircleAtPoint = bsxfun(@plus, SortValue(PlotFeatures(i))*[CircleX; CircleY], [a;b]);
    
    
    
%     if ismember( FeaturePoint3Dd3Dc.Location(SortOrder(i)), AllMatches2 )
        
        plot(CircleAtPoint(2,:), CircleAtPoint(1,:),'r')
%         
%     elseif ismember( FeaturePoint3Dd3Dc.Location(SortOrder(i)), AllMatches2 )
% %         
%         plot(CircleAtPoint(2,:), CircleAtPoint(1,:),'g')
% %         
%     else
% %         i
%         t=t+1;
% %         FeaturePoint2Dd3Dc.Location(SortOrder(i))
%         plot(CircleAtPoint(2,:), CircleAtPoint(1,:),'g')
%     end

    

end
t


NumPoints = length(FeaturePoint3Dd3Dc.Location);


AllMatches = zeros(NumPoints,1);
for i = 1 : NumPoints
    
    Matches = find( (FeaturePoint3Dd3Dc.Location(i) == FeaturePoint2Dd2Dc.Location) & (FeaturePoint3Dd3Dc.Scale(i) == FeaturePoint2Dd2Dc.Scale) );
    
    if isempty(Matches)
        i
    else
        if length(Matches)>1
            length(Matches)
        end
        AllMatches(i,1) = FeaturePoint3Dd3Dc.Location(i);
    end
    
end
AllMatches(AllMatches==0) = [];


AllMatches2 = zeros(NumPoints,1);
for i = 1 : NumPoints
    
    Matches = find( (FeaturePoint3Dd3Dc.Location(i) == FeaturePoint2Dd3Dc.Location) & (FeaturePoint3Dd3Dc.Scale(i) == FeaturePoint2Dd3Dc.Scale) );
    
    if isempty(Matches)
        i
    else
        if length(Matches)>1
            length(Matches)
        end
        AllMatches2(i,1) = FeaturePoint3Dd3Dc.Location(i);
    end
    
end
AllMatches2(AllMatches2==0) = [];





AllMatches2 = zeros(NumPoints,1);
for i = 1 : NumPoints
    
    
    [InterestPointA, InterestPointB] = ind2sub([MaxImageSize(1),MaxImageSize(2)], FeaturePoint3Dd3Dc.Location(i));
    
    [MatchSetA, MatchSetB] = ind2sub([MaxImageSize(1),MaxImageSize(2)], FeaturePoint2Dd2Dc.Location);
    
    Matches1 = find( sqrt( bsxfun(@minus, MatchSetA, InterestPointA).^2 + bsxfun(@minus, MatchSetB, InterestPointB).^2 ) < 1.1);
    

    
%     Matches1 = find( (FeaturePoint3Dd3Dc.Location(i) == FeaturePoint2Dd2Dc.Location) );
        
    Matches2 = find( abs(FeaturePoint3Dd3Dc.Scale(i) - FeaturePoint2Dd2Dc.Scale(Matches1))<10*eps );
    
    MatchesA = Matches1(Matches2);
    
    if length(MatchesA) > 1
        
        [MinVal, MinLoc] = min(sqrt( (MatchSetA(Matches1) -  InterestPointA).^2 + (MatchSetB(Matches1) -  InterestPointB).^2 ) );
        
        MatchesA = Matches1(MinLoc);
    end
    
    if isempty(MatchesA) && ~isempty(Matches1)
        MatchesScale = FeaturePoint2Dd2Dc.Scale(Matches1);
        MatchesB=[];
        for j = 1 : length(MatchesScale)
            
            ScaleVecMatch = find(ScaleParamterImage == MatchesScale(j));
            
            ScaleVecMatchNeighbor = [ScaleVecMatch-2:ScaleVecMatch+2];
            ScaleVecMatchNeighbor(ScaleVecMatchNeighbor < 1) = [];
            ScaleVecMatchNeighbor(ScaleVecMatchNeighbor > length(ScaleParamterImage)) = [];
            
            
            if any( FeaturePoint3Dd3Dc.Scale(i) == ScaleParamterImage(ScaleVecMatchNeighbor))
                MatchesB(j) = i;
            end
            
            
        end
        MatchesB(MatchesB ==0) = [];
        MatchesB = unique(MatchesB);
        
        if ~isempty(MatchesB)
            if length(MatchesB) > 1
                i
                length(MatchesB)
                break;
            end
            
            Matches = MatchesB;
        end
    else
        
        Matches = MatchesA;
        
    end
    
    
    if isempty(Matches)
%         i
    else
        if length(Matches)>1
            i
            length(Matches)
            break
        end
        AllMatches2(i,1) = FeaturePoint3Dd3Dc.Location(i);
    end
    
end
AllMatches2(AllMatches2==0) = [];
size(AllMatches2)





CP = find( FeaturePoint2Dd2Dc.Location == sub2ind([MaxImageSize(1),MaxImageSize(2)],312,153))


CS = FeaturePoint2Dd2Dc.Scale(CP)

[a,b] = ind2sub([MaxImageSize(1),MaxImageSize(2)], 1:NumPixels);

PointInCircle = find( sqrt( bsxfun(@minus, a, 312).^2 + bsxfun(@minus, b, 153).^2 ) <= CS)

IndPointInCircle = sub2ind([MaxImageSize(1),MaxImageSize(2)], a(PointInCircle), b(PointInCircle));


IntensityValues = abs(acosd(Image(IndPointInCircle))' - acosd(Image(FeaturePoint2Dd2Dc.Location(CP)))')


IntensityValues = acosd(Image(IndPointInCircle))'

IntensityValues = abs(IntensityValues - mean(IntensityValues));

p=hist(IntensityValues, 0:1:90)
p = p./max(p);
plot(p)








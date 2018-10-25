









function Features = buildHistogram(PointCloud, Keypoint)


if ~isfield(Keypoint, 'Count')
    Keypoint.Count = length(Keypoint.Location);
end


BinSize45 = 2 / 45;

BinSizeTheta = 180 / 45;

Tree = KDTreeSearcher(PointCloud.Location, 'Distance', 'Euclidean');

% Initialize variables
Features.HistAlpha = zeros(Keypoint.Count, 45);
Features.HistPhi = zeros(Keypoint.Count, 45);
Features.HistBeta = zeros(Keypoint.Count, 45);
Features.HistTheta = zeros(Keypoint.Count, 45);
Features.Scale = zeros(Keypoint.Count, 1);
Features.Location = zeros(Keypoint.Count, 3);
Features.LocationIndex = zeros(Keypoint.Count, 1);
Features.Normal = zeros(Keypoint.Count, 3);


warning('on','all')
WaitBar = waitbar(0, sprintf('Building Histogram for Keypoint %i of %i', 0, Keypoint.Count));

searcher = @(Tree, Point, Radius) rangesearch(Tree, Point, Radius);
Position = cellfun(searcher, repmat({Tree},Keypoint.Count,1), mat2cell(Keypoint.Location, ones(Keypoint.Count,1), 3), num2cell(Keypoint.Scale));


for i = 1 : Keypoint.Count

%     [Position, ~] = rangesearch(Tree, PointCloud.Location(NMSKeypoint.Location(i),:), NMSKeypoint.Scale(i));
    
    CurrentPosition = Position{i,1};
    CurrentPosition(1) = [];
    
    AllVertexDist = bsxfun(@minus, PointCloud.Location(CurrentPosition,:), Keypoint.Location(i,:));
    NormVertexDist = sqrt(sum(AllVertexDist.^2,2));
    
    %% Darboux Frame for each point
    u = Keypoint.Normal(i,:);
    v = bsxfun(@cross, bsxfun(@rdivide, AllVertexDist, NormVertexDist)', u')';
    v = bsxfun(@rdivide, v, sqrt(sum(v.^2,2)));
    w = bsxfun(@cross, u', v')';
    
    %% Angle Deviations
    % Beta and Phi are very similar to eachother. Will only need one.
    % 
    cosAlpha = sum( PointCloud.Normal(CurrentPosition,:) .* v, 2);
    cosPhi = bsxfun(@rdivide, AllVertexDist, NormVertexDist) * u';
    Theta = atan2d( sum(PointCloud.Normal(CurrentPosition,:) .* w, 2) , PointCloud.Normal(CurrentPosition,:) * u');
%     cosPsi = PointCloud.Normal(CurrentPosition,:) * u';
    cosBeta = sum( bsxfun(@rdivide, AllVertexDist, NormVertexDist) .* PointCloud.Normal(CurrentPosition,:), 2);
    sigma = Keypoint.Scale(i);
    
    
    %% Remove points from thin surfaces.
    % These points are on the other side of an object
%     RemoveIndex = cosPsi < 0.2; % Equivalent to 78.5 deg difference
%     if nnz(RemoveIndex)
%         warning('buildHistogram:: Removing points from thin surfaces')
%         disp(nnz(RemoveIndex))
%         disp(nnz(RemoveIndex)/length(RemoveIndex))
%         cosAlpha(RemoveIndex) = [];
%         cosPhi(RemoveIndex) = [];
%         Theta(RemoveIndex) = [];
% %         cosPsi(RemoveIndex) = [];
%         cosBeta(RemoveIndex) = [];
%     end
    
    %% Accumulate Angles in Bin
    % Need to implement a bilinear/trilinear interpolation scheme
    
    
    AlphaBins = histcounts(cosAlpha, (0:45)*BinSize45 - 1 );
    PhiBins = histcounts(cosPhi, (0:45)*BinSize45 - 1 );
    BetaBins = histcounts(cosBeta, (0:45)*BinSize45 - 1 );
    ThetaBins = histcounts(Theta, (0:45)*BinSizeTheta - 90 );

    
%     for j = 1 : 45
%         CurrentBin = j * BinSize45;
%         
%         % Alpha
%         AlphaInBin = find(cosAlpha >= (CurrentBin - BinSize45 - 1) & cosAlpha < (CurrentBin - 1));
%         AlphaBins(i,j) = size(AlphaInBin, 1);
%         
%         % Phi
%         PhiInBin = find(cosPhi >= (CurrentBin - BinSize45 - 1) & cosPhi < (CurrentBin - 1));
%         PhiBins(i,j) = size(PhiInBin, 1);
%         
%         % Psi
% %         PsiInBin = find(cosPsi >= (CurrentBin - BinSize45 - 1) & cosPsi < (CurrentBin - 1));
% %         PsiBins(i,j) = size(PsiInBin, 1);
%         
%         % Beta
%         BetaInBin = find(cosBeta >= (CurrentBin - BinSize45 - 1) & cosBeta < (CurrentBin - 1));
%         BetaBins(i,j) = size(BetaInBin, 1);
%         
%         
%         % Theta
%         CurrentThetaBin = j * BinSizeTheta;
%         ThetaInBin = find(Theta >= (CurrentThetaBin - BinSizeTheta - 90) & Theta < (CurrentThetaBin - 90));
%         ThetaBins(i,j) = size(ThetaInBin, 1);
%         
%     end
    
    if nnz(AlphaBins)
        AlphaBins = AlphaBins ./ sum(AlphaBins);
    end
    if nnz(PhiBins)
       PhiBins= PhiBins ./ sum(PhiBins); 
    end
%     if nnz(PsiBins(i,:))
%         PsiBins(i,:) = PsiBins(i,:) ./ sum(PsiBins(i,:));
%     end
    if nnz(BetaBins)
       BetaBins = BetaBins ./ sum(BetaBins); 
    end
    if nnz(ThetaBins)
        ThetaBins = ThetaBins ./ sum(ThetaBins);    
    end
    
    
    Features.HistAlpha(i,:) = AlphaBins;
    Features.HistPhi(i,:) = PhiBins;
%     Features.HistPsi(i,:) = PsiBins(i,:);
    Features.HistBeta(i,:) = BetaBins;
    Features.HistTheta(i,:) = ThetaBins;
    Features.Scale(i,1) = sigma;
    Features.Location(i,1:3) = Keypoint.Location(i,:);
    Features.LocationIndex(i,1) = Keypoint.LocationIndex(i,1);
    Features.Normal(i,1:3) = Keypoint.Normal(i,:);
    
%     figure
%     plot(1:45*4, [AlphaBins(i,:),PhiBins(i,:),BetaBins(i,:),ThetaBins(i,:) ] ,'-*')
%     drawnow
%     pause


waitbar(i/Keypoint.Count, WaitBar, sprintf('Building Histogram for Keypoint %i of %i', i, Keypoint.Count));

end

waitbar(i/Keypoint.Count, WaitBar, sprintf('Building Histograms Complete'));
close(WaitBar)

Features.Count = Keypoint.Count;

Features.Resolution = PointCloud.Resolution;

% 
% plot(1:45*5, [AlphaBins,PhiBins,PsiBins,BetaBins,ThetaBins ] ,'-*')
% plot(1:45, AlphaBins,'-*')
% plot(1:45, PhiBins,'-*')
% plot(1:45, PsiBins,'-*')
% plot(1:45, BetaBins,'-*')
% plot(1:45, ThetaBins,'-*')
% 








end


















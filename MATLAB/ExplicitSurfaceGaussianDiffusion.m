





KDTree = KDTreeSearcher(PointCloud.Location);

AverageEdgeLength = findAverageEdgeLengths(PointCloud.Location, PointCloud.Face, 1);



PointCloud.VertexArea = accumarray( reshape(PointCloud.Face, 3*PointCloud.FaceCount ,1), repmat(PointCloud.FaceArea,3,1), [PointCloud.LocationCount, 1] )/3;



count = 1;
maxcount = 0;
WaitBar = waitbar(0, sprintf('Vertex %i of %i', 0, PointCloud.LocationCount-1));
for i = 1 : PointCloud.LocationCount
    
     
    [Neigh, Dist] = rangesearch(KDTree, PointCloud.Location(i,:),2*sqrt(2)*AverageEdgeLength);
     
    CurrentNeighbors = Neigh{:};
    
    NumNeighbors = length(CurrentNeighbors);
    
    maxcount = maxcount + NumNeighbors;
    
    Position1(count:maxcount) = i * ones(NumNeighbors,1);
    Position2(count:maxcount) = CurrentNeighbors;
    
    Weights(count:maxcount) = 1/(2*pi*AverageEdgeLength^2) * PointCloud.VertexArea(i) * PointCloud.VertexArea(CurrentNeighbors) .* exp(- sum(bsxfun(@minus,PointCloud.Location(i,:), PointCloud.Location(CurrentNeighbors,:)).^2,2) / (2*AverageEdgeLength^2)  );


    count = maxcount + 1;
    waitbar(i/PointCloud.LocationCount, WaitBar, sprintf('Vertex %i of %i', i, PointCloud.LocationCount-1));
    
end

waitbar(i/PointCloud.LocationCount, WaitBar, sprintf('Building Complete'));
close(WaitBar)



Weights(Position1==0) = [];
Position2(Position1==0) = [];
Position1(Position1==0) = [];

MaxImageSize = [30 30];

NumVertex = prod(MaxImageSize);



A = sparse(Position1, Position2, Weights, NumVertex, NumVertex);

B = sparse(1:NumVertex,1:NumVertex,1./PointCloud.VertexArea) * A;

G = B ./ max(sum(B,2));


SignalExplicit = zeros(PointCloud.LocationCount, NumSteps);
SignalExplicit(:,1) = PointCloud.Signal;

NumSteps = 50;

WaitBar = waitbar(0, sprintf('Vertex %i of %i', 0, NumSteps-1));
figure
for i = 1 : NumSteps - 1
    
    SignalExplicit(:,i+1) = G * SignalExplicit(:,i);
    
    if flag
        flag
    end
    
    
     imshow(reshape(SignalExplicit(:,i), MaxImageSize(1), MaxImageSize(2)),[])
     drawnow
     
    waitbar(i/NumSteps, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumSteps-1));
    
end
waitbar(i/PointCloud.LocationCount, WaitBar, sprintf('Building Complete'));
close(WaitBar)





RadialDist3D = sqrt(sum(bsxfun(@minus, PointCloud.Location, PointCloud.Location(PointCloudCenter,:)).^2,2));

xgauss = 0:0.01:20;
figure
for i = 1 : NumSteps
    
    
    clf
    plot(RadialDist3D, SignalExplicit(:,i),'r.')
    hold on
    axis([0 10 0 2*max(SignalExplicit(:,i))])
    
    t = (i-1) * AverageEdgeLength^2;
    sqrt(t)
    
    Gauss = exp(-(xgauss.^2) / (2*t));
    Gauss = Gauss * max(SignalExplicit(:,i));
    plot(xgauss, Gauss,'b-')
    
    pause
end










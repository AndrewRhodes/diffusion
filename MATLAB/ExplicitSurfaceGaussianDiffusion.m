




AverageEdgeLength = findAverageEdgeLengths(PointCloud.Location, PointCloud.Face, 1);

 G = makeMeshGaussian(PointCloud, AverageEdgeLength, 2*sqrt(2)*AverageEdgeLength, 1);


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










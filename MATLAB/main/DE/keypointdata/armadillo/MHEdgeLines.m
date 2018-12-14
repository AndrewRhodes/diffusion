

% 0.1e-6, 0.25e-6, and 0.5e-6


thresh = 2.5e-6; % 10e-5; % 0.5e-6; % 0.1e-6; % 0.25e-6; % 
k=350;



figure
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), Signal(:,k), 'EdgeColor', 'none');
% view(180, -90)
view(0, 90)
ax = gca;
axis equal
axis off
hold on
cbar = colorbar('FontSize',15);


% CurrentDoGMax = DoGMaximum.Location(DoGMaximum.Level==k,1);
CurrentDoGMax = find(AbsDoG(:,k) < thresh);

Keep = true(length(CurrentDoGMax),1);

for i = 1 : length(CurrentDoGMax)
    
   CurrentNeighbors = Neighbors.Distance{CurrentDoGMax(i),1};
   CurrentValue = DoG(CurrentDoGMax(i),k);
   
   SurroundingValuesCurrent = AbsDoG(CurrentNeighbors,k); 
   
   AllSign = sign([SurroundingValuesCurrent', CurrentValue]);
   
   SignProd = movprod(AllSign,2,2);
   
   if ~( any(SignProd < 0) && any(SignProd > 0) )
        Keep(i) = 0;
   end
       
end

NewDoGMax = CurrentDoGMax(Keep);


plot3(PointCloud.Location(NewDoGMax,1), PointCloud.Location(NewDoGMax,2),...
    PointCloud.Location(NewDoGMax,3), 'k.','MarkerSize', 10)

plot3(PointCloud.Location(CurrentDoGMax,1), PointCloud.Location(CurrentDoGMax,2),...
    PointCloud.Location(CurrentDoGMax,3), 'k.','MarkerSize', 10)

ylabel(cbar, 'Mean Curvature', 'FontSize',30);
% ax.XLim = [-240 240];
% ax.YLim = [-115 200];
ax.XLim = [-90 70]; % bunny
ax.YLim = [15 135]; % bunny
ax.ZLim = [-50 50]; % bunny
% title(sprintf('Diffusion Step: %d, Absolute Scale: %0.2f, AbsDoG thesh: %0.1e',k, ScaleParameterAbsolute(k+1), thresh),'fontsize',25)
title(sprintf('Diffusion Step: %d, AbsDoG thesh: %0.1e',k, thresh),'fontsize',25)
caxis([-0.15,0.2]) % Armadillo
caxis([-0.1202, 0.1793]) % bunny
























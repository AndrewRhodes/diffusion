


[maxval, maxind] = max(DoG,[],2);

[minval, minind] = min(DoG,[],2);

ZeroMax = (maxind <= 87) | (maxind == 2999);

ZeroMin = (maxind <= 87) | (maxind == 2999);

LocationsMin = 1:PointCloud.LocationCount;
LocationsMax = 1:PointCloud.LocationCount;
LocationsMin(ZeroMin) = [];
LocationsMax(ZeroMax) = [];


MaxIndRemove = maxind(~ZeroMax);
MaxValRemove = maxval(~ZeroMax);

MinIndRemove = minind(~ZeroMin);
MinValRemove = minval(~ZeroMin);


MinQuant = quantile(MinValRemove, [0.25,0.5,0.75]);
MaxQuant = quantile(MaxValRemove, [0.25,0.5,0.75]);


% MinKeep = MinValRemove < (MinQuant(3) - 1.5*(MinQuant(3) - MinQuant(1)));
% MaxKeep = MaxValRemove > (MaxQuant(3) + 1.5*(MaxQuant(3) - MaxQuant(1)));

MinKeep = MinValRemove > MinQuant(3);
MaxKeep = MaxValRemove < MaxQuant(1);
% 
MinKeep = (MinValRemove < MinQuant(3)) & (MinValRemove > MinQuant(1));
MaxKeep = (MaxValRemove < MaxQuant(3)) & (MaxValRemove > MaxQuant(1));

LocationsMin = LocationsMin(MinKeep);
MinIndRemove = MinIndRemove(MinKeep);
MinValRemove = MinValRemove(MinKeep);

LocationsMax = LocationsMax(MaxKeep);
MaxIndRemove = MaxIndRemove(MaxKeep);
MaxValRemove = MaxValRemove(MaxKeep);


figure
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), Signal(:,1), 'EdgeColor', 'none');
view(180, -90)
ax = gca;
axis equal
axis off
hold on
cbar = colorbar('FontSize',15);
   
plot3(PointCloud.Location(LocationsMin,1), PointCloud.Location(LocationsMin,2),...
    PointCloud.Location(LocationsMin,3), 'k.','MarkerSize', 12)

plot3(PointCloud.Location(LocationsMax,1), PointCloud.Location(LocationsMax,2),...
    PointCloud.Location(LocationsMax,3), 'r.','MarkerSize', 12)

ylabel(cbar, 'Mean Curvature', 'FontSize', 30);
ax.XLim = [-240 240];
ax.YLim = [-115 200];
caxis([-0.15,0.2])

figure
boxplot(MaxIndRemove)

figure
boxplot(MinValRemove)
figure
boxplot(MaxValRemove)




% [a,b]=min(sqrt(sum((PointCloud.Location(Keypoint.Location,:) - [-11.29, 47.28, -11.14]).^2,2)))

kp = 4991;
% kp = 4139;
kp2 = 3682; %3813;

figure
% plot(ScaleParameter(2:end),DoG(Keypoint.Location(kp),:),'b.','linewidth', 4)
% hold on
% plot(ScaleParameter(2:end),DoG(Keypoint.Location(kp2),:),'r.','linewidth', 4)
% plot(Keypoint.Scale(kp)-PointCloud.ResolutionLocal(Keypoint.Location(kp)), DoG(Keypoint.Location(kp),Keypoint.Level(kp)),'b.','markersize',50)
% plot(Keypoint.Scale(kp2)-PointCloud.ResolutionLocal(Keypoint.Location(kp2)), DoG(Keypoint.Location(kp2),Keypoint.Level(kp2)),'r.','markersize',50)

plot(ScaleParameter(2:end)+PointCloud.ResolutionLocal(Keypoint.Location(kp)),DoG(Keypoint.Location(kp),:),'b.','linewidth', 4)
hold on
plot(ScaleParameter(2:end)+PointCloud.ResolutionLocal(Keypoint.Location(kp2)),DoG(Keypoint.Location(kp2),:),'r.','linewidth', 4)
plot(Keypoint.Scale(kp), DoG(Keypoint.Location(kp),Keypoint.Level(kp)),'b.','markersize',50)
plot(Keypoint.Scale(kp2), DoG(Keypoint.Location(kp2),Keypoint.Level(kp2)),'r.','markersize',50)
xy = [ScaleParameter(2)+min(PointCloud.ResolutionLocal(Keypoint.Location([kp,kp2]))), ScaleParameter(end)+max(PointCloud.ResolutionLocal(Keypoint.Location([kp,kp2])))];

plot([xy(1), xy(2)],[0,0],'k--' ,'linewidth', 2)
xlim([xy(1), xy(2)])
xlabel('Levels')
ylabel('DoG Response')
ax = gca;
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;


[SphereX, SphereY, SphereZ] = sphere(100);
SphereX = reshape(SphereX,[],1);
SphereY = reshape(SphereY,[],1);
SphereZ = reshape(SphereZ,[],1);


fig = figure('Units', 'Normalized', 'OuterPosition', [0 0 0.98 0.98]);
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), Signal(:,1), 'EdgeColor', 'none');
view(180, -90)
ax = gca;
axis equal
axis off
hold on

plot3(PointCloud.Location(Keypoint.Location,1), PointCloud.Location(Keypoint.Location,2),...
          PointCloud.Location(Keypoint.Location,3), 'k.','MarkerSize', 10)

SphereAtPoint = bsxfun(@plus, Keypoint.Scale(kp)*[SphereX, SphereY, SphereZ], PointCloud.Location(Keypoint.Location(kp),:));
hh1 = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
set(hh1, 'FaceColor','m', 'EdgeColor','none', 'FaceAlpha',0.8)

SphereAtPoint = bsxfun(@plus, Keypoint.Scale(kp2)*[SphereX, SphereY, SphereZ], PointCloud.Location(Keypoint.Location(kp2),:));
hh2 = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
set(hh2, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.8)

caxis([-0.15,0.2])

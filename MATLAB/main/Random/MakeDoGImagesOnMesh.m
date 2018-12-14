





global ProjectRoot;

ModelFolder = 'armadillo/';
Model = 'Armadillo_e1_100000';

FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat(ModelFolder,Model,'.ply');
FileNameModelOff = strcat(ModelFolder,Model,'.off');

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );

% PointCloud.Location = PointCloud.Location.*1.9868;

PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
PointCloud = findMeshResolution(PointCloud, 'Model');
PointCloud = findMeshNormals(PointCloud);


%%%%%%%%%%%%%%%%%%%%
% Load some signal too
%%%%%%%%%%%%%%%%%%%%%%%%


k1 = 50;
k2 = 100;




figure
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), Signal(:,k1), 'EdgeColor', 'none');
% view(180, -90)
view(0, 90)
ax = gca;
axis equal
axis off
hold on
% cbar = colorbar('FontSize',15);
ax.XLim = [-90 70];
ax.YLim = [15 135];
ax.ZLim = [-50 50];
caxis([-0.1202, 0.1793])



figure
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), Signal(:,k2), 'EdgeColor', 'none');
% view(180, -90)
view(0, 90)
ax = gca;
axis equal
axis off
hold on
% cbar = colorbar('FontSize',15);
ax.XLim = [-90 70];
ax.YLim = [15 135];
ax.ZLim = [-50 50];
caxis([-0.1202, 0.1793])

% 
% level = 2961; 
% %303;
% 330,844,2961
% 
% CurrentZeropoints = Zeropoint.Location( Zeropoint.Level == level);
% p1 = plot3(PointCloud.Location(CurrentZeropoints,1), PointCloud.Location(CurrentZeropoints,2),...
%     PointCloud.Location(CurrentZeropoints,3), 'k.','MarkerSize', 10);
% 
% 
% 
% 
% CurrentNewKeypoints = NewKeypoint.Location(NewKeypoint.Level == level);
% CurrentNewScales = NewKeypoint.Scale(NewKeypoint.Level == level);
% CurrentSigns = NewKeypoint.Sign(NewKeypoint.Level == level);
% 
% for j = 1 : length(CurrentNewKeypoints)
%     SphereAtPoint = bsxfun(@plus, CurrentNewScales(j)*[SphereX, SphereY, SphereZ], PointCloud.Location(CurrentNewKeypoints(j),:));
%     hh = surfl(reshape(SphereAtPoint(:,1),101,101), reshape(SphereAtPoint(:,2),101,101), reshape(SphereAtPoint(:,3),101,101));
%     set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)
%     
%     if CurrentSigns(j) < 0
%         set(hh, 'FaceColor','b', 'EdgeColor','none', 'FaceAlpha',0.5)
%     else
%         set(hh, 'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.5)
%     end
% end
% drawnow
% pause

    
    
    

SigDif = (Signal(:,k2) - Signal(:,k1)) * (ScaleParameter(k1,1)^2 / ( ScaleParameter(k2,1)^2 - ScaleParameter(k1,1)^2 ) );

Quant = quantile(SigDif, [0.1,0.9]);


figure
trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), SigDif, 'EdgeColor', 'none');
% view(180, -90)
view(0, 90)
ax = gca;
axis equal
axis off
hold on
% cbar = colorbar('FontSize',15);
ax.XLim = [-90 70];
ax.YLim = [15 135];
ax.ZLim = [-50 50];
caxis([Quant(1), Quant(2)])


























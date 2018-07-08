

close all
clear
clc

global ProjectRoot; % Additional Paths

NumIter = 50;
Use4D = 0;
UseNMS = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Keypoint Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelData = {'bunny','armadillo','buddha','dragon'};
ModelNames = {'Bunny_e1', 'Armadillo_e1', 'Buddha_e1_50000', 'Dragon_e1_50000'};
modelDataLength = length(modelData);
Error = cell(modelDataLength,modelDataLength);
Match = cell(modelDataLength,modelDataLength);
NoMatch = cell(modelDataLength,modelDataLength);
MultipleMatch = cell(modelDataLength,modelDataLength);

for k = 1 : modelDataLength
    k
    Model = strcat(modelData{k},'/', ModelNames{k});
    FileLocationModel = strcat(ProjectRoot,'/models/object/');
    FileNameModelPly = strcat(Model,'.ply');
    
    [PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
    
    PointCloud.LocationCount = size(PointCloud.Location,1);
    PointCloud.FaceCount = size(PointCloud.Face, 1);
    PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
    PointCloud = findMeshResolution(PointCloud, 'Model');
%     PointCloud = findMeshNormals(PointCloud);
%     NormalRotations = findNormalsRotation(PointCloud.Normal);
    


    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',modelData{k},'/Std_0.1/');
    [Error{k,1}, Match{k,1}, NoMatch{k,1}, MultipleMatch{k,1}] = findKeypointError(PointCloud, FileLocation, NumIter, UseNMS, Use4D);
    
    
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',modelData{k},'/Std_0.2/');
    [Error{k,2}, Match{k,2}, NoMatch{k,2}, MultipleMatch{k,2}] = findKeypointError(PointCloud, FileLocation, NumIter, UseNMS, Use4D);
    
    
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',modelData{k},'/Std_0.3/');
    [Error{k,3}, Match{k,3}, NoMatch{k,3}, MultipleMatch{k,3}] = findKeypointError(PointCloud, FileLocation, NumIter, UseNMS, Use4D);
    
    
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',modelData{k},'/Std_0.4/');
    [Error{k,4}, Match{k,4}, NoMatch{k,4}, MultipleMatch{k,4}] = findKeypointError(PointCloud, FileLocation, NumIter, UseNMS, Use4D);
    
    clear PointCloud
    
end



figure(1)
bp1 = boxplot([Error{1,1}.ScaleRepeat; Error{2,1}.ScaleRepeat; Error{3,1}.ScaleRepeat; Error{4,1}.ScaleRepeat],...
    [ones(NumIter,1); 1.25*ones(NumIter,1); 1.5*ones(NumIter,1); 1.75*ones(NumIter,1)],'positions',[1,1.25,1.5,1.75]);
hold on
bp2 = boxplot([Error{1,2}.ScaleRepeat; Error{2,2}.ScaleRepeat; Error{3,2}.ScaleRepeat; Error{4,2}.ScaleRepeat],...
    [3*ones(NumIter,1); 3.25*ones(NumIter,1); 3.5*ones(NumIter,1); 3.75*ones(NumIter,1)],'positions',[3,3.25,3.5,3.75]);
bp3 = boxplot([Error{1,3}.ScaleRepeat; Error{2,3}.ScaleRepeat; Error{3,3}.ScaleRepeat; Error{4,3}.ScaleRepeat],...
    [5*ones(NumIter,1); 5.25*ones(NumIter,1); 5.5*ones(NumIter,1); 5.75*ones(NumIter,1)],'positions',[5,5.25,5.5,5.75]);
bp4 = boxplot([Error{1,4}.ScaleRepeat; Error{2,4}.ScaleRepeat; Error{3,4}.ScaleRepeat; Error{4,4}.ScaleRepeat],...
    [7*ones(NumIter,1); 7.25*ones(NumIter,1); 7.5*ones(NumIter,1); 7.75*ones(NumIter,1)],'positions',[7,7.25,7.5,7.75]);
axis([0,8,0,1])
set(bp1,'Linewidth',3,'Markersize',15)
set(bp2,'Linewidth',3,'Markersize',15)
set(bp3,'Linewidth',3,'Markersize',15)
set(bp4,'Linewidth',3,'Markersize',15)
ylim([0 1])
xlabel('Noise')
ylabel('Repeatability')
ax = gca;
ax.XTick = [1.375, 3.375, 5.375, 7.375];
ax.XTickLabel = {'0.1','0.2','0.3','0.4'};
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;


figure(2)
bp1 = boxplot([Error{1,1}.RelativeRepeat; Error{2,1}.RelativeRepeat; Error{3,1}.RelativeRepeat; Error{4,1}.RelativeRepeat],...
    [ones(NumIter,1); 1.25*ones(NumIter,1); 1.5*ones(NumIter,1); 1.75*ones(NumIter,1)],'positions',[1,1.25,1.5,1.75]);
hold on
bp2 = boxplot([Error{1,2}.RelativeRepeat; Error{2,2}.RelativeRepeat; Error{3,2}.RelativeRepeat; Error{4,2}.RelativeRepeat],...
    [3*ones(NumIter,1); 3.25*ones(NumIter,1); 3.5*ones(NumIter,1); 3.75*ones(NumIter,1)],'positions',[3,3.25,3.5,3.75]);
bp3 = boxplot([Error{1,3}.RelativeRepeat; Error{2,3}.RelativeRepeat; Error{3,3}.RelativeRepeat; Error{4,3}.RelativeRepeat],...
    [5*ones(NumIter,1); 5.25*ones(NumIter,1); 5.5*ones(NumIter,1); 5.75*ones(NumIter,1)],'positions',[5,5.25,5.5,5.75]);
bp4 = boxplot([Error{1,4}.RelativeRepeat; Error{2,4}.RelativeRepeat; Error{3,4}.RelativeRepeat; Error{4,4}.RelativeRepeat],...
    [7*ones(NumIter,1); 7.25*ones(NumIter,1); 7.5*ones(NumIter,1); 7.75*ones(NumIter,1)],'positions',[7,7.25,7.5,7.75]);
axis([0,8,0,1])
set(bp1,'Linewidth',3,'Markersize',15)
set(bp2,'Linewidth',3,'Markersize',15)
set(bp3,'Linewidth',3,'Markersize',15)
set(bp4,'Linewidth',3,'Markersize',15)
ylim([0 1])
xlabel('Noise')
ylabel('Repeatability')
ax = gca;
ax.XTick = [1.375, 3.375, 5.375, 7.375];
ax.XTickLabel = {'0.1','0.2','0.3','0.4'};
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;



figure(3)
bp1 = boxplot([Error{1,1}.Count; Error{2,1}.Count; Error{3,1}.Count; Error{4,1}.Count],...
    [ones(NumIter,1); 1.25*ones(NumIter,1); 1.5*ones(NumIter,1); 1.75*ones(NumIter,1)],'positions',[1,1.25,1.5,1.75]);
hold on
bp2 = boxplot([Error{1,2}.Count; Error{2,2}.Count; Error{3,2}.Count; Error{4,2}.RelativeRepeat],...
    [3*ones(NumIter,1); 3.25*ones(NumIter,1); 3.5*ones(NumIter,1); 3.75*ones(NumIter,1)],'positions',[3,3.25,3.5,3.75]);
bp3 = boxplot([Error{1,3}.Count; Error{2,3}.Count; Error{3,3}.Count; Error{4,3}.RelativeRepeat],...
    [5*ones(NumIter,1); 5.25*ones(NumIter,1); 5.5*ones(NumIter,1); 5.75*ones(NumIter,1)],'positions',[5,5.25,5.5,5.75]);
bp4 = boxplot([Error{1,4}.Count; Error{2,4}.Count; Error{3,4}.Count; Error{4,4}.RelativeRepeat],...
    [7*ones(NumIter,1); 7.25*ones(NumIter,1); 7.5*ones(NumIter,1); 7.75*ones(NumIter,1)],'positions',[7,7.25,7.5,7.75]);
set(bp1,'Linewidth',3,'Markersize',15)
set(bp2,'Linewidth',3,'Markersize',15)
set(bp3,'Linewidth',3,'Markersize',15)
set(bp4,'Linewidth',3,'Markersize',15)
ylims = get(gca, 'YLim');
set(gca, 'Ylim', [0 ylims(2)]);
xlabel('Noise')
ax = gca;
ax.XTick = [1.375, 3.375, 5.375, 7.375];
ax.XTickLabel = {'0.1','0.2','0.3','0.4'};
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;




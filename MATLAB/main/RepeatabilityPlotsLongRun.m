

close all
clear
clc

global ProjectRoot; % Additional Paths

NumIter = 1;
Use4D = 0;
UseNMS = 1;
t_scale = 0.7;

tauFraction = 1/10;
NumSteps = 250;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Keypoint Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelData = {'bunny','armadillo','buddha','dragon'};
ModelNames = {'Bunny_e1', 'Armadillo_e1_100000', 'Buddha_e1_50000', 'Dragon_e1_50000'};
modelDataLength = length(modelData);
Error = cell(modelDataLength,modelDataLength);
Match = cell(modelDataLength,modelDataLength);
NoMatch = cell(modelDataLength,modelDataLength);
MultipleMatch = cell(modelDataLength,modelDataLength);
NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];
NoiseVecLength = length(NoiseVec);


for k = 1 : modelDataLength
    k
    Model = strcat(modelData{k},'/', ModelNames{k})
    FileLocationModel = strcat(ProjectRoot,'/models/object/');
    FileNameModelPly = strcat(Model,'.ply');
    
    [PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
    
    PointCloud.LocationCount = size(PointCloud.Location,1);
    PointCloud.FaceCount = size(PointCloud.Face, 1);
    PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
    PointCloud = findMeshResolution(PointCloud, 'Model');
    PointCloud = findMeshNormals(PointCloud);
    
    scale_min = sqrt(2*NumSteps*PointCloud.Resolution*tauFraction);
    
    for j = 1 : NoiseVecLength

        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',modelData{k},'/LongRun/Std_',num2str(NoiseVec(j)),'/');
        [Error{k,j}, Match{k,j}, NoMatch{k,j}, MultipleMatch{k,j}] = findKeypointError(PointCloud, FileLocation, NumIter, UseNMS, Use4D, t_scale, scale_min);

    end
    
    clear PointCloud
    
end


positions = [1, 1.25, 1.5, 1.75];

figure(1)
hold on
for j = 1 : NoiseVecLength
    cp = positions + (j-1)*2;
    bp = boxplot([Error{1,j}.ScaleRepeat; Error{2,j}.ScaleRepeat; Error{3,j}.ScaleRepeat; Error{4,j}.ScaleRepeat],...
    [cp(1)*ones(NumIter,1); cp(2)*ones(NumIter,1); cp(3)*ones(NumIter,1); cp(4)*ones(NumIter,1)],'positions',cp);
    set(bp,'Linewidth',2,'Markersize',8)
end

% bp1 = boxplot([Error{1,1}.ScaleRepeat; Error{2,1}.ScaleRepeat; Error{3,1}.ScaleRepeat; Error{4,1}.ScaleRepeat],...
%     [ones(NumIter,1); 1.25*ones(NumIter,1); 1.5*ones(NumIter,1); 1.75*ones(NumIter,1)],'positions',[1,1.25,1.5,1.75]);
% hold on
% bp2 = boxplot([Error{1,2}.ScaleRepeat; Error{2,2}.ScaleRepeat; Error{3,2}.ScaleRepeat; Error{4,2}.ScaleRepeat],...
%     [3*ones(NumIter,1); 3.25*ones(NumIter,1); 3.5*ones(NumIter,1); 3.75*ones(NumIter,1)],'positions',[3,3.25,3.5,3.75]);
% bp3 = boxplot([Error{1,3}.ScaleRepeat; Error{2,3}.ScaleRepeat; Error{3,3}.ScaleRepeat; Error{4,3}.ScaleRepeat],...
%     [5*ones(NumIter,1); 5.25*ones(NumIter,1); 5.5*ones(NumIter,1); 5.75*ones(NumIter,1)],'positions',[5,5.25,5.5,5.75]);
% bp4 = boxplot([Error{1,4}.ScaleRepeat; Error{2,4}.ScaleRepeat; Error{3,4}.ScaleRepeat; Error{4,4}.ScaleRepeat],...
%     [7*ones(NumIter,1); 7.25*ones(NumIter,1); 7.5*ones(NumIter,1); 7.75*ones(NumIter,1)],'positions',[7,7.25,7.5,7.75]);
axis([0,8,0,1])
% set(bp1,'Linewidth',2,'Markersize',8)
% set(bp2,'Linewidth',2,'Markersize',8)
% set(bp3,'Linewidth',2,'Markersize',8)
% set(bp4,'Linewidth',2,'Markersize',8)
ylim([0 1])
xlabel('Noise')
ylabel('Percentage')
ax = gca;
ax.XTick = mean(positions) + ((1:NoiseVecLength)-1)*2;
ax.XTickLabel = num2cell(NoiseVec);%{'0.1','0.2','0.3','0.4', '0.5'};
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;


figure(2)
hold on
for j = 1 : NoiseVecLength
    cp = positions + (j-1)*2;
    bp = boxplot([Error{1,j}.RelativeRepeat; Error{2,j}.RelativeRepeat; Error{3,j}.RelativeRepeat; Error{4,j}.RelativeRepeat],...
    [cp(1)*ones(NumIter,1); cp(2)*ones(NumIter,1); cp(3)*ones(NumIter,1); cp(4)*ones(NumIter,1)],'positions',cp);
    set(bp,'Linewidth',2,'Markersize',8)
end
% bp1 = boxplot([Error{1,1}.RelativeRepeat; Error{2,1}.RelativeRepeat; Error{3,1}.RelativeRepeat; Error{4,1}.RelativeRepeat],...
%     [ones(NumIter,1); 1.25*ones(NumIter,1); 1.5*ones(NumIter,1); 1.75*ones(NumIter,1)],'positions',[1,1.25,1.5,1.75]);
% hold on
% bp2 = boxplot([Error{1,2}.RelativeRepeat; Error{2,2}.RelativeRepeat; Error{3,2}.RelativeRepeat; Error{4,2}.RelativeRepeat],...
%     [3*ones(NumIter,1); 3.25*ones(NumIter,1); 3.5*ones(NumIter,1); 3.75*ones(NumIter,1)],'positions',[3,3.25,3.5,3.75]);
% bp3 = boxplot([Error{1,3}.RelativeRepeat; Error{2,3}.RelativeRepeat; Error{3,3}.RelativeRepeat; Error{4,3}.RelativeRepeat],...
%     [5*ones(NumIter,1); 5.25*ones(NumIter,1); 5.5*ones(NumIter,1); 5.75*ones(NumIter,1)],'positions',[5,5.25,5.5,5.75]);
% bp4 = boxplot([Error{1,4}.RelativeRepeat; Error{2,4}.RelativeRepeat; Error{3,4}.RelativeRepeat; Error{4,4}.RelativeRepeat],...
%     [7*ones(NumIter,1); 7.25*ones(NumIter,1); 7.5*ones(NumIter,1); 7.75*ones(NumIter,1)],'positions',[7,7.25,7.5,7.75]);
axis([0,8,0,1])
% set(bp1,'Linewidth',2,'Markersize',8)
% set(bp2,'Linewidth',2,'Markersize',8)
% set(bp3,'Linewidth',2,'Markersize',8)
% set(bp4,'Linewidth',2,'Markersize',8)
ylim([0 1])
xlabel('Noise')
ylabel('Percentage')
ax = gca;
ax.XTick = mean(positions) + ((1:NoiseVecLength)-1)*2;
ax.XTickLabel = {'0.1','0.2','0.3','0.4', '0.5'};
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;



figure(3)
hold on
for j = 1 : NoiseVecLength
    cp = positions + (j-1)*2;
    bp = boxplot([Error{1,j}.Count; Error{2,j}.Count; Error{3,j}.Count; Error{4,j}.Count],...
    [cp(1)*ones(NumIter,1); cp(2)*ones(NumIter,1); cp(3)*ones(NumIter,1); cp(4)*ones(NumIter,1)],'positions',cp);
    set(bp,'Linewidth',2,'Markersize',8)
end
% bp1 = boxplot([Error{1,1}.Count; Error{2,1}.Count; Error{3,1}.Count; Error{4,1}.Count],...
%     [ones(NumIter,1); 1.25*ones(NumIter,1); 1.5*ones(NumIter,1); 1.75*ones(NumIter,1)],'positions',[1,1.25,1.5,1.75]);
% hold on
% bp2 = boxplot([Error{1,2}.Count; Error{2,2}.Count; Error{3,2}.Count; Error{4,2}.RelativeRepeat],...
%     [3*ones(NumIter,1); 3.25*ones(NumIter,1); 3.5*ones(NumIter,1); 3.75*ones(NumIter,1)],'positions',[3,3.25,3.5,3.75]);
% bp3 = boxplot([Error{1,3}.Count; Error{2,3}.Count; Error{3,3}.Count; Error{4,3}.RelativeRepeat],...
%     [5*ones(NumIter,1); 5.25*ones(NumIter,1); 5.5*ones(NumIter,1); 5.75*ones(NumIter,1)],'positions',[5,5.25,5.5,5.75]);
% bp4 = boxplot([Error{1,4}.Count; Error{2,4}.Count; Error{3,4}.Count; Error{4,4}.RelativeRepeat],...
%     [7*ones(NumIter,1); 7.25*ones(NumIter,1); 7.5*ones(NumIter,1); 7.75*ones(NumIter,1)],'positions',[7,7.25,7.5,7.75]);
% set(bp1,'Linewidth',2,'Markersize',8)
% set(bp2,'Linewidth',2,'Markersize',8)
% set(bp3,'Linewidth',2,'Markersize',8)
% set(bp4,'Linewidth',2,'Markersize',8)
ylims = get(gca, 'YLim');
set(gca, 'Ylim', [0 ylims(2)]);
xlim([0 8])
ylabel('Number')
xlabel('Noise')
ax = gca;
ax.XTick = mean(positions) + ((1:NoiseVecLength)-1)*2;
ax.XTickLabel = {'0.1','0.2','0.3','0.4', '0.5'};
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;






close all
clear
clc

global ProjectRoot; % Additional Paths

NumIter = 1;
Use4D = 0;
UseNMS = 1;
t_scale = 0.7;
level_min = 250;
t_range = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Keypoint Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelData = {'bunny','armadillo','buddha','dragon','itokawa'}; % {'armadillo'}; %
ModelNames = {'Bunny_e1', 'Armadillo_e1_100000', 'Buddha_e1_50000', 'Dragon_e1_50000','Itokawa_e1_80000'}; % {'Armadillo_e1_100000'}; %
modelDataLength = length(modelData);
Error = cell(modelDataLength, modelDataLength);
Match = cell(modelDataLength, modelDataLength);
NoMatch = cell(modelDataLength, modelDataLength);
MultipleMatch = cell(modelDataLength, modelDataLength);
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
    %     PointCloud = findMeshNormals(PointCloud);
    %     NormalRotations = findNormalsRotation(PointCloud.Normal);
    %     [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);
    
    load(strcat(ModelNames{k},'_Neighbors.mat'),'Neighbors')
    PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);
    
    for kk = 1 : NoiseVecLength
        
        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',modelData{k},'/LongRun/Std_',num2str(NoiseVec(kk)),'/');
        [Error{k,kk}, NoMatch{k,kk}] = findKeypointErrorNew(PointCloud, FileLocation, NumIter, UseNMS, t_scale, level_min, t_range);
        %         [Error{k,kk}, Match{k,kk}, NoMatch{k,kk}, MultipleMatch{k,kk}] = findKeypointError(PointCloud, FileLocation, NumIter, UseNMS, Use4D, t_scale, level_min);
        
    end
    
    clear PointCloud
    
end


positions = [1, 1.2, 1.4, 1.6, 1.8];
% positions = [1, 1.25, 1.5, 1.75];

if NumIter ~= 1
    figure(1)
    hold on
    for kk = 1 : NoiseVecLength
        cp = positions + (kk-1)*2;
        bp = boxplot([Error{1,kk}.ScaleRepeat; Error{2,kk}.ScaleRepeat; Error{3,kk}.ScaleRepeat; Error{4,kk}.ScaleRepeat; Error{5,kk}.ScaleRepeat],...
            [cp(1)*ones(NumIter,1); cp(2)*ones(NumIter,1); cp(3)*ones(NumIter,1); cp(4)*ones(NumIter,1); cp(5)*ones(NumIter,1)],'positions',cp);
        set(bp,'Linewidth',2,'Markersize',8)
    end
    axis([0,8,0,1])
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
    for kk = 1 : NoiseVecLength
        cp = positions + (kk-1)*2;
        bp = boxplot([Error{1,kk}.RelativeRepeat; Error{2,kk}.RelativeRepeat; Error{3,kk}.RelativeRepeat; Error{4,kk}.RelativeRepeat; Error{5,kk}.RelativeRepeat],...
            [cp(1)*ones(NumIter,1); cp(2)*ones(NumIter,1); cp(3)*ones(NumIter,1); cp(4)*ones(NumIter,1); cp(5)*ones(NumIter,1)],'positions',cp);
        set(bp,'Linewidth',2,'Markersize',8)
    end
    axis([0,8,0,1])
    
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
    for kk = 1 : NoiseVecLength
        cp = positions + (kk-1)*2;
        bp = boxplot([Error{1,kk}.Count; Error{2,kk}.Count; Error{3,kk}.Count; Error{4,kk}.Count; Error{5,kk}.Count],...
            [cp(1)*ones(NumIter,1); cp(2)*ones(NumIter,1); cp(3)*ones(NumIter,1); cp(4)*ones(NumIter,1); cp(5)*ones(NumIter,1)],'positions',cp);
        set(bp,'Linewidth',2,'Markersize',8)
    end
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
    
    
else
    
    
    figure(1)
    hold on
    y=[];
    for kk = 1 : NoiseVecLength
        y = [y; Error{1,kk}.ScaleRepeat, Error{2,kk}.ScaleRepeat, Error{3,kk}.ScaleRepeat, Error{4,kk}.ScaleRepeat, Error{5,kk}.ScaleRepeat];
    end
    bar(NoiseVec, y)
    xlabel('Noise')
    ylabel('Percentage')
    ax = gca;
    ax.XTick = NoiseVec;
    ax.YTick = [0,0.25,0.5,0.75,1];
    ax.XAxis.FontSize = 50;
    ax.YAxis.FontSize = 50;
    ylim([0 1])
    
    
    figure(2)
    hold on
    y=[];
    for kk = 1 : NoiseVecLength
        y = [y; Error{1,kk}.RelativeRepeat, Error{2,kk}.RelativeRepeat, Error{3,kk}.RelativeRepeat, Error{4,kk}.RelativeRepeat, Error{5,kk}.RelativeRepeat];
    end
    bar(NoiseVec,y)
    xlabel('Noise')
    ylabel('Percentage')
    ax = gca;
    ax.XTick = NoiseVec;
    ax.YTick = [0,0.25,0.5,0.75,1];
    ax.XAxis.FontSize = 50;
    ax.YAxis.FontSize = 50;
    ylim([0 1])
    
    figure(3)
    hold on
    y=[];
    for kk = 1 : NoiseVecLength
        y = [y; Error{1,kk}.Count, Error{2,kk}.Count, Error{3,kk}.Count, Error{4,kk}.Count, Error{5,kk}.Count];
    end
    bar(NoiseVec, y)
    ax = gca;
    ax.XTick = NoiseVec;
    ax.XAxis.FontSize = 50;
    ax.YAxis.FontSize = 50;
%     ylims = get(gca, 'YLim');
%     set(ax, 'Ylim', [0 ylims(2)]);
    ylabel('Number')
    xlabel('Noise')
    
end
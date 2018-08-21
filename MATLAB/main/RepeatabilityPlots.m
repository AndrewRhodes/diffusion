

close all
clear
clc

global ProjectRoot; % Additional Paths

NumIter = 50;
Use4D = 0;
UseNMS = 1;
t_scale = 0.7;
level_min = 0;
t_range = 3;
NoiseType = 'Signal'; % 'Signal', 'Vertex'
lines = 1;

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
    
    if NumIter == 1
        Location = strcat('/main/DE/keypointdata/',modelData{k},'/',NoiseType,'Noise/LongRun/Std_');
    elseif NumIter == 50
        Location = strcat('/main/DE/keypointdata/',modelData{k},'/',NoiseType,'Noise/ShortRun/Std_');
    end
    
    for kk = 1 : NoiseVecLength
        
        FileLocation = strcat(ProjectRoot,Location,num2str(NoiseVec(kk)),'/');
        [Error{k,kk}, NoMatch{k,kk}] = findKeypointErrorNew(PointCloud, FileLocation, NumIter, UseNMS, t_scale, level_min, t_range);
        %         [Error{k,kk}, Match{k,kk}, NoMatch{k,kk}, MultipleMatch{k,kk}] = findKeypointError(PointCloud, FileLocation, NumIter, UseNMS, Use4D, t_scale, level_min);
        
    end
    
    clear PointCloud
    
end

positions = [1, 1.2, 1.4, 1.6, 1.8];
% positions = [1, 1.25, 1.5, 1.75];

if NumIter == 50
    
    figure(1)
    hold on
    for kk = 1 : NoiseVecLength
        cp = positions + (kk-1)*2;
        bp = boxplot([Error{1,kk}.ScaleRepeat; Error{2,kk}.ScaleRepeat; Error{3,kk}.ScaleRepeat; Error{4,kk}.ScaleRepeat; Error{5,kk}.ScaleRepeat],...
            [cp(1)*ones(NumIter,1); cp(2)*ones(NumIter,1); cp(3)*ones(NumIter,1); cp(4)*ones(NumIter,1); cp(5)*ones(NumIter,1)],'positions',cp);
        set(bp,'Linewidth',2,'Markersize',8)
    end
    axis([0,10,0,1])
    ylim([0 1])
    xlabel('Noise')
    ylabel('Percentage')
    ax = gca;
    ax.YTick = [0,0.2,0.4,0.6,0.8,1];
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
    axis([0,10,0,1])
    
    ylim([0 1])
    xlabel('Noise')
    ylabel('Percentage')
    ax = gca;
    ax.YTick = [0,0.2,0.4,0.6,0.8,1];
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
    xlim([0 10])
    ylabel('Number')
    xlabel('Noise')
    ax = gca;
    %     ax.YTick = [0,1000,2000,3000,4000,5000,6000];
    ax.XTick = mean(positions) + ((1:NoiseVecLength)-1)*2;
    ax.XTickLabel = {'0.1','0.2','0.3','0.4', '0.5'};
    ax.XAxis.FontSize = 50;
    ax.YAxis.FontSize = 50;
    
    
elseif NumIter == 1
    
    if lines
        markers = {'+-','o-','d-','s-','x-'};
        
        figure(1)
        hold on
        for kk = 1 : modelDataLength
            y = [Error{kk,1}.ScaleRepeat, Error{kk,2}.ScaleRepeat, Error{kk,3}.ScaleRepeat, Error{kk,4}.ScaleRepeat, Error{kk,5}.ScaleRepeat];
            plot(NoiseVec, y, markers{kk}, 'Linewidth', 7, 'MarkerSize', 20)
        end
        xlabel('Noise')
        ylabel('Percentage')
        ax = gca;
        ax.XTick = NoiseVec;
        ax.YTick = [0,0.2,0.4,0.6,0.8,1];
        ax.XAxis.FontSize = 50;
        ax.YAxis.FontSize = 50;
        ylim([0 1])
        legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'FontSize',45,'Orientation','Horizontal')
        
        
        
        figure(2)
        hold on
        for kk = 1 : modelDataLength
            y = [Error{kk,1}.RelativeRepeat, Error{kk,2}.RelativeRepeat, Error{kk,3}.RelativeRepeat, Error{kk,4}.RelativeRepeat, Error{kk,5}.RelativeRepeat];
            plot(NoiseVec, y, markers{kk}, 'Linewidth', 7, 'MarkerSize', 20)
        end
        xlabel('Noise')
        ylabel('Percentage')
        ax = gca;
        ax.XTick = NoiseVec;
        ax.YTick = [0,0.2,0.4,0.6,0.8,1];
        ax.XAxis.FontSize = 50;
        ax.YAxis.FontSize = 50;
        ylim([0 1])
        legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'FontSize',45,'Orientation','Horizontal')
        
        
        
        figure(3)
        hold on
        for kk = 1 : modelDataLength
            y = [Error{kk,1}.Count, Error{kk,2}.Count, Error{kk,3}.Count, Error{kk,4}.Count, Error{kk,5}.Count];
            plot(NoiseVec, y, markers{kk}, 'Linewidth', 7, 'MarkerSize', 20)
        end
        xlabel('Noise')
        ylabel('Percentage')
        ax = gca;
        ax.XTick = NoiseVec;
        ax.XAxis.FontSize = 50;
        ax.YAxis.FontSize = 50;
        legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'FontSize',45,'Orientation','Horizontal')
        
        
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
        ax.YTick = [0,0.2,0.4,0.6,0.8,1];
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
        ax.YTick = [0,0.2,0.4,0.6,0.8,1];
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
    
    
end



KeypointLevelCounts = zeros(3000,modelDataLength);
% KeypointScaleCounts = cell(modelDataLength,1);

for k = 1 : modelDataLength
    
    if NumIter == 1
        Location = strcat('/main/DE/keypointdata/',modelData{k},'/',NoiseType,'Noise/LongRun/Std_');
    elseif NumIter == 50
        Location = strcat('/main/DE/keypointdata/',modelData{k},'/',NoiseType,'Noise/ShortRun/Std_');
    end
    
     FileLocation = strcat(ProjectRoot,Location,'0.1','/');
     
    FileName = strcat('NMSKeypoint','.mat');
    load(fullfile(FileLocation, FileName),'NMSKeypoint')

    LevelTable = tabulate(NMSKeypoint.Level);
    
    KeypointLevelCounts(LevelTable(:,1),k) = LevelTable(:,2);
    
%     ScaleTable = tabulate(NMSKeypoint.Scale);
%     KeypointScaleCounts{k,1} = ScaleTable(:,1:2);
    
end




markers = {'+-','o-','d-','s-','x-'};

figure(4)
hold on
for kk = 1 : modelDataLength
    plot(1:3000, KeypointLevelCounts(:,kk), '.-', 'Linewidth', 7, 'MarkerSize', 20)
end
xlim([0,150])
xlabel('Level')
ylabel('Count')
ax = gca;
ax.XAxis.FontSize = 50;
ax.YAxis.FontSize = 50;













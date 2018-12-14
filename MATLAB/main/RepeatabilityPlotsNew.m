

close all
clear
clc

global ProjectRoot; % Additional Paths

NumIter = 10;
t_scale = 0.7;
level_min = 2;
t_range = 3/4;
NoiseType = 'Signal'; % 'Signal', 'Vertex'
lines = 1;
ByModel = 1;

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
    load(strcat(ModelNames{k},'_Neighbors.mat'), 'Neighbors')
    PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);

    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',modelData{k},'/',NoiseType,'Noise/NewRun/');
    
    for kk = 1 : NoiseVecLength
        
        NoiseFileLocation = strcat(FileLocation,'Std_',num2str(NoiseVec(kk)),'/');
        [Error{k,kk}, NoMatch{k,kk}] = findKeypointErrorNewer(PointCloud, FileLocation, NoiseFileLocation, NumIter, t_scale, level_min, t_range);

    end
    
    clear PointCloud
    
end

positions = [1, 1.2, 1.4, 1.6, 1.8];
% positions = [1, 1.25, 1.5, 1.75];

if NumIter > 1
    
    
    
    figure
    hold on
    for kk = 1 : NoiseVecLength
        cp = positions + (kk-1)*2;
        Data = [];
        Points = [];
        for jj = 1 : modelDataLength
            Data = [Data; Error{jj,kk}.ScaleRepeat];
            Points = [Points; cp(jj)*ones(NumIter,1)];
        end
        bp = boxplot(Data, Points, 'positions',cp(1:modelDataLength));
        set(bp,'Linewidth',2,'Markersize',8)        
%         bp = boxplot([Error{1,kk}.ScaleRepeat; Error{2,kk}.ScaleRepeat; Error{3,kk}.ScaleRepeat; Error{4,kk}.ScaleRepeat; Error{5,kk}.ScaleRepeat],...
%             [cp(1)*ones(NumIter,1); cp(2)*ones(NumIter,1); cp(3)*ones(NumIter,1); cp(4)*ones(NumIter,1); cp(5)*ones(NumIter,1)],'positions',cp);
    end
    axis([0,2*NoiseVecLength,0,1])
    ylim([0 1])
    xlabel('Noise')
    ylabel('Percentage')
    ax = gca;
    ax.YTick = [0,0.2,0.4,0.6,0.8,1];
    ax.XTick = mean(positions) + ((1:NoiseVecLength)-1)*2;
    ax.XTickLabel = num2cell(NoiseVec);%{'0.1','0.2','0.3','0.4', '0.5'};
    ax.XAxis.FontSize = 50;
    ax.YAxis.FontSize = 50;
    
    
    
    figure
    hold on
    for kk = 1 : NoiseVecLength
        cp = positions + (kk-1)*2;
        Data = [];
        Points = [];
        for jj = 1 : modelDataLength
            Data = [Data; Error{jj,kk}.RelativeRepeat];
            Points = [Points; cp(jj)*ones(NumIter,1)];
        end
        bp = boxplot(Data, Points, 'positions',cp(1:modelDataLength));
        set(bp,'Linewidth',2,'Markersize',8)
%         bp = boxplot([Error{1,kk}.RelativeRepeat; Error{2,kk}.RelativeRepeat; Error{3,kk}.RelativeRepeat; Error{4,kk}.RelativeRepeat; Error{5,kk}.RelativeRepeat],...
%             [cp(1)*ones(NumIter,1); cp(2)*ones(NumIter,1); cp(3)*ones(NumIter,1); cp(4)*ones(NumIter,1); cp(5)*ones(NumIter,1)],'positions',cp);
    end
    axis([0,2*NoiseVecLength,0,1])
    
    ylim([0 1])
    xlabel('Noise')
    ylabel('Percentage')
    ax = gca;
    ax.YTick = [0,0.2,0.4,0.6,0.8,1];
    ax.XTick = mean(positions) + ((1:NoiseVecLength)-1)*2;
    ax.XTickLabel = {'0.1','0.2','0.3','0.4', '0.5'};
    ax.XAxis.FontSize = 50;
    ax.YAxis.FontSize = 50;
    
    
    
    figure
    hold on
    for kk = 1 : NoiseVecLength
        cp = positions + (kk-1)*2;
        Data = [];
        Points = [];
        for jj = 1 : modelDataLength
            Data = [Data; Error{jj,kk}.Count];
            Points = [Points; cp(jj)*ones(NumIter,1)];
        end
        bp = boxplot(Data, Points, 'positions',cp(1:modelDataLength));
        set(bp,'Linewidth',2,'Markersize',8)
%         bp = boxplot([Error{1,kk}.Count; Error{2,kk}.Count; Error{3,kk}.Count; Error{4,kk}.Count; Error{5,kk}.Count],...
%             [cp(1)*ones(NumIter,1); cp(2)*ones(NumIter,1); cp(3)*ones(NumIter,1); cp(4)*ones(NumIter,1); cp(5)*ones(NumIter,1)],'positions',cp);
        ylims(kk,:) = get(gca, 'YLim');
    end
        
    set(gca, 'Ylim', [0 max(ylims(:,2))]);
    xlim([0 2*NoiseVecLength])
    ylabel('Count')
    xlabel('Noise')
    ax = gca;
    %     ax.YTick = [0,1000,2000,3000,4000,5000,6000];
    ax.XTick = mean(positions) + ((1:NoiseVecLength)-1)*2;
    ax.XTickLabel = {'0.1','0.2','0.3','0.4', '0.5'};
    ax.YTick = [0,40,80,120,160];
    ax.XAxis.FontSize = 50;
    ax.YAxis.FontSize = 50;
    
    
    if ByModel
        
        figure
        hold on
        for jj = 1 : modelDataLength
            cp = positions + (jj -1)*2;
            Data = [];
            Points = [];
            for kk = 1 : NoiseVecLength
                Data = [Data; Error{jj,kk}.ScaleRepeat];
                Points = [Points; cp(kk)*ones(NumIter,1)];
            end
            bp = boxplot(Data, Points, 'positions', cp);
            set(bp,'Linewidth',2,'Markersize',8)
        end
        
        axis([0,2*modelDataLength,0,1])
        ylim([0 1])
        %     xlabel('Model')
        ylabel('Percentage')
        ax = gca;
        ax.YTick = [0,0.2,0.4,0.6,0.8,1];
        ax.XTick = 1.25 + ((1:modelDataLength)-1)*2;
        ax.XTickLabel = modelData;
        ax.XAxis.FontSize = 50;
        ax.YAxis.FontSize = 50;
        xtickangle(ax,-30)
        title('Scale Repeat','FontSize',50)
        
        
        figure
        hold on
        for jj = 1 : modelDataLength
            cp = positions + (jj -1)*2;
            Data = [];
            Points = [];
            for kk = 1 : NoiseVecLength
                Data = [Data; Error{jj,kk}.RelativeRepeat];
                Points = [Points; cp(kk)*ones(NumIter,1)];
            end
            bp = boxplot(Data, Points, 'positions', cp);
            set(bp,'Linewidth',2,'Markersize',8)
        end
        
        axis([0,2*modelDataLength,0,1])
        ylim([0 1])
        %     xlabel('Model')
        ylabel('Percentage')
        ax = gca;
        ax.YTick = [0,0.2,0.4,0.6,0.8,1];
        ax.XTick = 1.25 + ((1:modelDataLength)-1)*2;
        ax.XTickLabel = modelData;
        ax.XAxis.FontSize = 50;
        ax.YAxis.FontSize = 50;
        xtickangle(ax,-30)
        title('Relative Repeat','FontSize',50)
        
        
        
        figure
        hold on
        for jj = 1 : modelDataLength
            cp = positions + (jj -1)*2;
            Data = [];
            Points = [];
            for kk = 1 : NoiseVecLength
                Data = [Data; Error{jj,kk}.Count];
                Points = [Points; cp(kk)*ones(NumIter,1)];
            end
            bp = boxplot(Data, Points, 'positions', cp);
            set(bp,'Linewidth',2,'Markersize',8)
            ylims(jj,:) = get(gca, 'YLim');
        end
        
        set(gca, 'Ylim', [0 max(ylims(:,2))]);
        xlim([0,2*modelDataLength])
        ylabel('Count')
        ax = gca;
        ax.YTick = [0,40,80,120,160];
        ax.XTick = 1.25 + ((1:modelDataLength)-1)*2;
        ax.XTickLabel = modelData;
        ax.XAxis.FontSize = 50;
        ax.YAxis.FontSize = 50;
        xtickangle(ax,-30)
        title('Absolute Repeat','FontSize',50)
       
    end
    
    
    
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
        ylabel('Count')
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
        ylabel('Count')
        xlabel('Noise')
        
        
        
    end
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analyzing different step sizes
KeypointLevelCounts = zeros(3000,modelDataLength);
FileLocationModel = strcat(ProjectRoot,'/models/object/');
FileNameModelPly = strcat('itokawa/Itokawa_e1_80000','.ply');

[PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
PointCloud = findMeshResolution(PointCloud, 'Model');
PointCloud

    
figure(4)
hold on
for k = 1 : 3
    
    clear NMSKeypoint
    if k == 1
        Location = strcat('/main/DE/keypointdata/itokawa/SignalNoise/LongRun/Std_');
        FileLocation = strcat(ProjectRoot,Location,'0.1','/');
        FileName = strcat('NMSKeypoint','.mat');
        load(fullfile(FileLocation, FileName),'NMSKeypoint')
        c=10;
    elseif k == 2
        load('NMSKeypointEbar5.mat')
        c=5;
    elseif k == 3
        load('NMSKeypointEbar1.mat')
        c=1;
    end
    
    tau = PointCloud.Resolution / c;   
    ScaleParameter = findScaleParamter(tau, 1, 3000, 'Laplacian', 'Natural');
    
%     if k == 1 
%         load('NMSKeypointEbar20.mat')
%         c=20;
%     elseif k == 2
%         load('NMSKeypointEbar15.mat')
%         c=15;
%     elseif k == 3
%         Location = strcat('/main/DE/keypointdata/itokawa/SignalNoise/LongRun/Std_');
%         FileLocation = strcat(ProjectRoot,Location,'0.1','/');
%         FileName = strcat('NMSKeypoint','.mat');
%         load(fullfile(FileLocation, FileName),'NMSKeypoint')
%         c=10;
%     elseif k == 4
%         load('NMSKeypointEbar5.mat')
%         c=5;
%     elseif k == 5
%        load('NMSKeypointEbar1.mat')
%         c=1;
%     end
    
    LevelTable = tabulate(NMSKeypoint.Level);
    
    KeypointLevelCounts(LevelTable(:,1),k) = LevelTable(:,2);
    
    KeypointLevels(1:3000,k) = ScaleParameter ;
    
    plot(KeypointLevels(:,k), KeypointLevelCounts(:,k), '.-')
end
legend({'$\tau = \bar{e}/10$','$\tau = \bar{e}/5$','$\tau = \bar{e}/1$'},'interpreter','Latex','FontSize',30)
% legend({'c = 20','c = 15','c = 10','c = 5','c = 1'},'FontSize',30)
% legend({'$t=\bar{e}/10$','$t=\bar{e}/5$','$t=\bar{e}/1$'},'interpreter','Latex','FontSize',30)
xlim([0,5])
xlabel({'$\sigma_i$'},'interpreter','Latex')
ylabel('Count')
ax = gca;
ax.XAxis.FontSize = 25;
ax.YAxis.FontSize = 25;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


KeypointLevelCounts = zeros(3000,modelDataLength);
% KeypointScaleCounts = cell(modelDataLength,1);
markers = {'+-','o-','d-','s-','x-'};
a = [10,5,1];
figure(5)
hold on

for k = 1 : modelDataLength
    
    
    Model = strcat(modelData{k},'/', ModelNames{k})
    FileLocationModel = strcat(ProjectRoot,'/models/object/');
    FileNameModelPly = strcat(Model,'.ply');
    
    [PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
    PointCloud = findMeshResolution(PointCloud, 'Model');
    PointCloud
    
    tau = PointCloud.Resolution / 10;   
    
    ScaleParameter = findScaleParamter(tau, 1, 3000, 'Laplacian', 'Natural');

	Location = strcat('/main/DE/keypointdata/',modelData{k},'/SignalNoise/LongRun/Std_');

    
	FileLocation = strcat(ProjectRoot,Location,'0.1','/');
     
    FileName = strcat('NMSKeypoint','.mat');
    load(fullfile(FileLocation, FileName),'NMSKeypoint')

    for j = 1 : 
    LevelTable = tabulate(NMSKeypoint.Level);
    
    KeypointLevelCounts(LevelTable(:,1),k) = LevelTable(:,2);
    
    KeypointLevels(1:3000,k) = (ScaleParameter) ./ PointCloud.Resolution;

    plot(KeypointLevels(:,k), KeypointLevelCounts(:,k), '.-')
    
end
legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'Orientation','Horizontal') 
xlim([0, 5])
xlabel({'$n \tau / \bar{e}$'},'interpreter','Latex')
xlabel({'$\sigma_i / \bar{e}$'},'interpreter','Latex')
ylabel('Count')
ax = gca;
ax.XAxis.FontSize = 25;
ax.YAxis.FontSize = 25;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    
    

KeypointLevelCounts = zeros(3000,modelDataLength);
% KeypointScaleCounts = cell(modelDataLength,1);
markers = {'+-','o-','d-','s-','x-'};
a = [10,5,1];
figure(5)
hold on

for k = 1 : modelDataLength
    
    
    Model = strcat(modelData{k},'/', ModelNames{k})
    FileLocationModel = strcat(ProjectRoot,'/models/object/');
    FileNameModelPly = strcat(Model,'.ply');
    
    [PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
    PointCloud = findMeshResolution(PointCloud, 'Model');
    PointCloud
    
    tau = PointCloud.Resolution / 10;   
    
    ScaleParameter = findScaleParamter(tau, 1, 3000, 'Laplacian', 'Natural');

	Location = strcat('/main/DE/keypointdata/',modelData{k},'/SignalNoise/LongRun/Std_');

    
	FileLocation = strcat(ProjectRoot,Location,'0.1','/');
     
    FileName = strcat('NMSKeypoint','.mat');
    load(fullfile(FileLocation, FileName),'NMSKeypoint')

    
    LevelTable = tabulate(NMSKeypoint.Level);
    
    KeypointLevelCounts(LevelTable(:,1),k) = LevelTable(:,2);
    
    KeypointLevels(1:3000,k) = (ScaleParameter) ./ PointCloud.Resolution;

    plot(KeypointLevels(:,k), KeypointLevelCounts(:,k), '.-')
    
end
legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'Orientation','Horizontal') 
xlim([0, 5])
xlabel({'$n \tau / \bar{e}$'},'interpreter','Latex')
xlabel({'$\sigma_i / \bar{e}$'},'interpreter','Latex')
ylabel('Count')
ax = gca;
ax.XAxis.FontSize = 25;
ax.YAxis.FontSize = 25;



for kk = 1 : 3 %modelDataLength
    plot(KeypointLevels(:,kk)*PointCloud.Resolution / (PointCloud.Resolution/a(kk)), KeypointLevelCounts(:,kk), '.-')%, 'Linewidth', 7, 'MarkerSize', 20)
%     plot(1:3000, KeypointLevelCounts(:,kk), '.-')%, 'Linewidth', 7, 'MarkerSize', 20)
end
% legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'Orientation','Horizontal') 

xlim([0,150])
xlabel({'$\frac{Level}{\bar{e}}$'},'interpreter','Latex')
ylabel('Count')
ax = gca;
ax.XAxis.FontSize = 25;
ax.YAxis.FontSize = 25;
% legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'FontSize',45,'Orientation','Horizontal')        




% Find percentage of model size for a given number of steps
n = 100;
tauFraction = 1/10;
for k = 1 : modelDataLength
    
    Model = strcat(modelData{k},'/', ModelNames{k})
    FileLocationModel = strcat(ProjectRoot,'/models/object/');
    FileNameModelPly = strcat(Model,'.ply');
    
    [PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
    PointCloud = findMeshResolution(PointCloud, 'Model');
    PointCloud.Resolution
    
%     100*sqrt(2*n*PointCloud.Resolution*tauFraction)./range(PointCloud.Location)
end







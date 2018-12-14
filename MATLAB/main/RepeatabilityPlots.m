

close all
clear
clc

global ProjectRoot; % Additional Paths

UseNMS = 1;
t_scale = 1/sqrt(2);
level_min = 1;
t_range = 1/2;
NoiseType = 'Vertex'; % 'Signal', 'Vertex' %
RunType = 'Special'; % 'Umbrella', 'Cotangent', 'Mesh', 'Mesh_euc_te8', 'Mesh_te8' %
lines = 1;
ByModel = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Keypoint Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelData = {'buddha'}; %{'bunny','armadillo','buddha','dragon','itokawa'}; %{'dragon'}; % {'bunny','armadillo','buddha','dragon','itokawa'}; % {'buddha','dragon'}; %
ModelNames = {'Buddha_e1_50000'}; %{'Bunny_e1','Armadillo_e1_100000','Buddha_e1_50000', 'Dragon_e1_50000', 'Itokawa_e1_80000'}; %{'Dragon_e1_50000'}; %   {'Bunny_e1','Armadillo_e1_100000','Buddha_e1_50000', 'Dragon_e1_50000', 'Itokawa_e1_80000'};  % {'Buddha_e1_50000','Dragon_e1_50000'}; %

modelDataLength = length(modelData);
Error = cell(modelDataLength, modelDataLength);
Match = cell(modelDataLength, modelDataLength);
NoMatch = cell(modelDataLength, modelDataLength);
MultipleMatch = cell(modelDataLength, modelDataLength);
NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];
NoiseVecLength = length(NoiseVec);

% 
% ScaleParameter = sqrt((options.hs/2)^2*(1:NumSteps)');
% ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);
% [MinVal, MinLoc]=min(abs(ScaleParameter - 1 * PointCloud.Resolution))
% 
% 
% 
% Keypoint = KeypointB;
% ZeroOut = Keypoint.Level < MinLoc;
% Keypoint.Location(ZeroOut) = [];
% Keypoint.Level(ZeroOut) = [];
% Keypoint.Scale(ZeroOut) = [];
% Keypoint.Sign(ZeroOut) = [];
% Keypoint




for k = 1 : modelDataLength

    Model = strcat(modelData{k},'/', ModelNames{k})
    FileLocationModel = strcat(ProjectRoot,'/models/object/');
    FileNameModelPly = strcat(Model,'.ply');
    FileLocationNeighbors = strcat('/media/andrew/WDRhodes/diffusiondata/',modelData{k},'/neighbors/');
    
    
    [PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
    
    PointCloud.LocationCount = size(PointCloud.Location,1);
    PointCloud.FaceCount = size(PointCloud.Face, 1);
    PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
    PointCloud = findMeshResolution(PointCloud, 'Model');
    %     PointCloud = findMeshNormals(PointCloud);
    %     NormalRotations = findNormalsRotation(PointCloud.Normal);
    %     [Neighbors, NeighborFaces, PointCloud] = findAdjacentNeighbors(PointCloud);
    
    load(strcat(FileLocationNeighbors, ModelNames{k},'_Neighbors.mat'), 'Neighbors')
    PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);
    PointCloud = findMeshNormals(PointCloud);
    
    
    
    FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',modelData{k},'/',NoiseType,'Noise/',RunType,'/');
    
%     if strcmpi(NoiseType,'Vertex')
%         PountCloudNoiseModelLocation = strcat(ProjectRoot, '/models/object/',modelData{k},'/VertexNoise/');
%     end

    
    for kk = 1 : NoiseVecLength
        
        NoiseFileLocation = strcat(FileLocation,'Std_',num2str(NoiseVec(kk)),'/');
        
        if UseNMS
            NumIter = numel(dir(strcat(NoiseFileLocation, 'NMS*.mat')));
        else
            NumIter = numel(dir(strcat(NoiseFileLocation, 'Keypoint*.mat')));
        end
        
        NumIter %= NumIter/2
       
%         if strcmpi(NoiseType,'Vertex')
%             
%             [Error{k,kk}, NoMatch{k,kk}] = findKeypointErrorVertexNoise(PointCloud, NoiseVec(kk), ...
%                 FileLocation, NoiseFileLocation, NumIter, UseNMS, t_scale, level_min, t_range);
%             
%         else

            [Error{k,kk}, NoMatch{k,kk}] = findKeypointErrorNew(PointCloud, FileLocation,...
                NoiseFileLocation, NumIter, UseNMS, t_scale, level_min, t_range);
            
%             NoisyPCLocation = strcat('/home/andrew/GitAndrew/diffusion/MATLAB/models/object/',modelData{k},'/VertexNoise/',ModelNames{k},'_sigma',num2str(NoiseVec(kk)));
%             NoisyPCNeighbors = strcat(FileLocationNeighbors, ModelNames{k},'_Neighbors_sigma',num2str(NoiseVec(kk)));
% %     
%             [Error{k,kk}, NoMatch{k,kk}] = keypointErrorTest(PointCloud, FileLocation,...
%                 NoiseFileLocation, NumIter, UseNMS, t_scale, level_min, t_range, NoisyPCLocation, NoisyPCNeighbors);

            
            % [Error{k,kk}, Match{k,kk}, NoMatch{k,kk}, MultipleMatch{k,kk}] = findKeypointError(PointCloud, FileLocation, NumIter, UseNMS, Use4D, t_scale, level_min);
        
%         end
        
    end
    
    clear PointCloud
    
end

positions = [1, 1.2, 1.4, 1.6, 1.8];
% positions = [1, 1.25, 1.5, 1.75];

% if NumIter > 1
%     
%     figure(1)
%     hold on
%     for kk = 1 : NoiseVecLength
%         cp = positions + (kk-1)*2;
%         bp = boxplot([Error{1,kk}.ScaleRepeat; Error{2,kk}.ScaleRepeat; Error{3,kk}.ScaleRepeat; Error{4,kk}.ScaleRepeat; Error{5,kk}.ScaleRepeat],...
%             [cp(1)*ones(NumIter,1); cp(2)*ones(NumIter,1); cp(3)*ones(NumIter,1); cp(4)*ones(NumIter,1); cp(5)*ones(NumIter,1)],'positions',cp);
%         set(bp,'Linewidth',2,'Markersize',8)
%     end
%     axis([0,10,0,1])
%     ylim([0 1])
%     xlabel('Noise')
%     ylabel('Percentage')
%     ax = gca;
%     ax.YTick = [0,0.2,0.4,0.6,0.8,1];
%     ax.XTick = mean(positions) + ((1:NoiseVecLength)-1)*2;
%     ax.XTickLabel = num2cell(NoiseVec);%{'0.1','0.2','0.3','0.4', '0.5'};
%     ax.XAxis.FontSize = 50;
%     ax.YAxis.FontSize = 50;
%     
%     
%     figure(2)
%     hold on
%     for kk = 1 : NoiseVecLength
%         cp = positions + (kk-1)*2;
%         bp = boxplot([Error{1,kk}.RelativeRepeat; Error{2,kk}.RelativeRepeat; Error{3,kk}.RelativeRepeat; Error{4,kk}.RelativeRepeat; Error{5,kk}.RelativeRepeat],...
%             [cp(1)*ones(NumIter,1); cp(2)*ones(NumIter,1); cp(3)*ones(NumIter,1); cp(4)*ones(NumIter,1); cp(5)*ones(NumIter,1)],'positions',cp);
%         set(bp,'Linewidth',2,'Markersize',8)
%     end
%     axis([0,10,0,1])
%     
%     ylim([0 1])
%     xlabel('Noise')
%     ylabel('Percentage')
%     ax = gca;
%     ax.YTick = [0,0.2,0.4,0.6,0.8,1];
%     ax.XTick = mean(positions) + ((1:NoiseVecLength)-1)*2;
%     ax.XTickLabel = {'0.1','0.2','0.3','0.4', '0.5'};
%     ax.XAxis.FontSize = 50;
%     ax.YAxis.FontSize = 50;
%     
%     
%     
%     figure(3)
%     hold on
%     for kk = 1 : NoiseVecLength
%         cp = positions + (kk-1)*2;
%         bp = boxplot([Error{1,kk}.Count; Error{2,kk}.Count; Error{3,kk}.Count; Error{4,kk}.Count; Error{5,kk}.Count],...
%             [cp(1)*ones(NumIter,1); cp(2)*ones(NumIter,1); cp(3)*ones(NumIter,1); cp(4)*ones(NumIter,1); cp(5)*ones(NumIter,1)],'positions',cp);
%         set(bp,'Linewidth',2,'Markersize',8)
%         ylims(kk,:) = get(gca, 'YLim');
%     end
%     set(gca, 'Ylim', [0 max(ylims(:,2))]);
%     xlim([0 10])
%     ylabel('Count')
%     xlabel('Noise')
%     ax = gca;
%     %     ax.YTick = [0,1000,2000,3000,4000,5000,6000];
%     ax.XTick = mean(positions) + ((1:NoiseVecLength)-1)*2;
%     ax.XTickLabel = {'0.1','0.2','0.3','0.4', '0.5'};
%     ax.XAxis.FontSize = 50;
%     ax.YAxis.FontSize = 50;
%     
%     if ByModel
        
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
        
        xlim([0,2*modelDataLength])
        set(gca, 'Ylim', [0 max(ylims(:,2))]);
        ylabel('Count')
        ax = gca;
%         ax.YTick = [0,40,80,120,160];
        ax.XTick = 1.25 + ((1:modelDataLength)-1)*2;
        ax.XTickLabel = modelData;
        ax.XAxis.FontSize = 50;
        ax.YAxis.FontSize = 50;
        xtickangle(ax,-30)
        
%     end
%     
%     
%     
% elseif NumIter == 1
%     
%     if lines
%         markers = {'+-','o-','d-','s-','x-'};
%         
%         figure(1)
%         hold on
%         for kk = 1 : modelDataLength
%             y = [Error{kk,1}.ScaleRepeat, Error{kk,2}.ScaleRepeat, Error{kk,3}.ScaleRepeat, Error{kk,4}.ScaleRepeat, Error{kk,5}.ScaleRepeat];
%             plot(NoiseVec, y, markers{kk}, 'Linewidth', 7, 'MarkerSize', 20)
%         end
%         xlabel('Noise')
%         ylabel('Percentage')
%         ax = gca;
%         ax.XTick = NoiseVec;
%         ax.YTick = [0,0.2,0.4,0.6,0.8,1];
%         ax.XAxis.FontSize = 50;
%         ax.YAxis.FontSize = 50;
%         ylim([0 1])
%         legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'FontSize',45,'Orientation','Horizontal')
%         
%         
%         
%         figure(2)
%         hold on
%         for kk = 1 : modelDataLength
%             y = [Error{kk,1}.RelativeRepeat, Error{kk,2}.RelativeRepeat, Error{kk,3}.RelativeRepeat, Error{kk,4}.RelativeRepeat, Error{kk,5}.RelativeRepeat];
%             plot(NoiseVec, y, markers{kk}, 'Linewidth', 7, 'MarkerSize', 20)
%         end
%         xlabel('Noise')
%         ylabel('Percentage')
%         ax = gca;
%         ax.XTick = NoiseVec;
%         ax.YTick = [0,0.2,0.4,0.6,0.8,1];
%         ax.XAxis.FontSize = 50;
%         ax.YAxis.FontSize = 50;
%         ylim([0 1])
%         legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'FontSize',45,'Orientation','Horizontal')
%         
%         
%         
%         figure(3)
%         hold on
%         for kk = 1 : modelDataLength
%             y = [Error{kk,1}.Count, Error{kk,2}.Count, Error{kk,3}.Count, Error{kk,4}.Count, Error{kk,5}.Count];
%             plot(NoiseVec, y, markers{kk}, 'Linewidth', 7, 'MarkerSize', 20)
%         end
%         xlabel('Noise')
%         ylabel('Count')
%         ax = gca;
%         ax.XTick = NoiseVec;
%         ax.XAxis.FontSize = 50;
%         ax.YAxis.FontSize = 50;
%         legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'FontSize',45,'Orientation','Horizontal')
%         
%         
%     else
%         
%         
%         figure(1)
%         hold on
%         y=[];
%         for kk = 1 : NoiseVecLength
%             y = [y; Error{1,kk}.ScaleRepeat, Error{2,kk}.ScaleRepeat, Error{3,kk}.ScaleRepeat, Error{4,kk}.ScaleRepeat, Error{5,kk}.ScaleRepeat];
%         end
%         bar(NoiseVec, y)
%         xlabel('Noise')
%         ylabel('Percentage')
%         ax = gca;
%         ax.XTick = NoiseVec;
%         ax.YTick = [0,0.2,0.4,0.6,0.8,1];
%         ax.XAxis.FontSize = 50;
%         ax.YAxis.FontSize = 50;
%         ylim([0 1])
%         legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'FontSize',45,'Orientation','Horizontal')
%         
%         
%         figure(2)
%         hold on
%         y=[];
%         for kk = 1 : NoiseVecLength
%             y = [y; Error{1,kk}.RelativeRepeat, Error{2,kk}.RelativeRepeat, Error{3,kk}.RelativeRepeat, Error{4,kk}.RelativeRepeat, Error{5,kk}.RelativeRepeat];
%         end
%         bar(NoiseVec,y)
%         xlabel('Noise')
%         ylabel('Percentage')
%         ax = gca;
%         ax.XTick = NoiseVec;
%         ax.YTick = [0,0.2,0.4,0.6,0.8,1];
%         ax.XAxis.FontSize = 50;
%         ax.YAxis.FontSize = 50;
%         ylim([0 1])
%         legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'FontSize',45,'Orientation','Horizontal')
%         
%         
%         
%         
%         figure(3)
%         hold on
%         y=[];
%         for kk = 1 : NoiseVecLength
%             y = [y; Error{1,kk}.Count, Error{2,kk}.Count, Error{3,kk}.Count, Error{4,kk}.Count, Error{5,kk}.Count];
%         end
%         bar(NoiseVec, y)
%         ax = gca;
%         ax.XTick = NoiseVec;
%         ax.XAxis.FontSize = 50;
%         ax.YAxis.FontSize = 50;
%         %     ylims = get(gca, 'YLim');
%         %     set(ax, 'Ylim', [0 ylims(2)]);
%         ylabel('Count')
%         xlabel('Noise')
%         legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'FontSize',45,'Orientation','Horizontal')
%         
%         
%     end
%     
%     
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Analyzing different step sizes
% % KeypointLevelCounts = zeros(3000,modelDataLength);
% FileLocationModel = strcat(ProjectRoot,'/models/object/');
% FileNameModelPly = strcat('itokawa/Itokawa_e1_80000','.ply');
% 
% [PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
% PointCloud = findMeshResolution(PointCloud, 'Model');
% PointCloud
% 
%     
% figure(4)
% hold on
% for k = 1 : 3
%     
%     clear NMSKeypoint
%     if k == 1
%         Location = strcat('/main/DE/keypointdata/itokawa/SignalNoise/LongRun/Std_');
%         FileLocation = strcat(ProjectRoot,Location,'0.1','/');
%         FileName = strcat('NMSKeypoint','.mat');
%         load(fullfile(FileLocation, FileName),'NMSKeypoint')
%         c=10;
%     elseif k == 2
%         load('NMSKeypointEbar5.mat')
%         c=5;
%     elseif k == 3
%         load('NMSKeypointEbar1.mat')
%         c=1;
%     end
%     
%     tau = PointCloud.Resolution / c;   
%     ScaleParameter = findScaleParamter(tau, 1, 3000, 'Laplacian', 'Natural');
%     
% %     if k == 1 
% %         load('NMSKeypointEbar20.mat')
% %         c=20;
% %     elseif k == 2
% %         load('NMSKeypointEbar15.mat')
% %         c=15;
% %     elseif k == 3
% %         Location = strcat('/main/DE/keypointdata/itokawa/SignalNoise/LongRun/Std_');
% %         FileLocation = strcat(ProjectRoot,Location,'0.1','/');
% %         FileName = strcat('NMSKeypoint','.mat');
% %         load(fullfile(FileLocation, FileName),'NMSKeypoint')
% %         c=10;
% %     elseif k == 4
% %         load('NMSKeypointEbar5.mat')
% %         c=5;
% %     elseif k == 5
% %        load('NMSKeypointEbar1.mat')
% %         c=1;
% %     end
%     
%     LevelTable = tabulate(NMSKeypoint.Level);
%     KeypointLevelCounts = zeros(3000,modelDataLength);
%     KeypointLevelCounts(LevelTable(:,1),k) = LevelTable(:,2);
%     
%     KeypointLevels(1:3000,k) = ScaleParameter ;
%     
%     plot(KeypointLevels(:,k), KeypointLevelCounts(:,k), '.-')
% end
% legend({'$\tau = \bar{e}/10$','$\tau = \bar{e}/5$','$\tau = \bar{e}/1$'},'interpreter','Latex','FontSize',30)
% % legend({'c = 20','c = 15','c = 10','c = 5','c = 1'},'FontSize',30)
% % legend({'$t=\bar{e}/10$','$t=\bar{e}/5$','$t=\bar{e}/1$'},'interpreter','Latex','FontSize',30)
% xlim([0,5])
% xlabel({'$\sigma_i$'},'interpreter','Latex')
% ylabel('Count')
% ax = gca;
% ax.XAxis.FontSize = 25;
% ax.YAxis.FontSize = 25;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% KeypointLevelCounts = zeros(3000,modelDataLength);
% % KeypointScaleCounts = cell(modelDataLength,1);
% % markers = {'+-','o-','d-','s-','x-'};
% markers = {'-',':','-.','--','-o'};
% figure
% hold on
% 
% for k = 1 : modelDataLength
%     
%     
%     Model = strcat(modelData{k},'/', ModelNames{k})
%     FileLocationModel = strcat(ProjectRoot,'/models/object/');
%     FileNameModelPly = strcat(Model,'.ply');
%     
%     [PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
%     PointCloud = findMeshResolution(PointCloud, 'Model');
%     PointCloud
%     
%     tau = PointCloud.Resolution / 10;   
%     
%     ScaleParameter = findScaleParamter(tau, 1, 3000, 'Laplacian', 'Natural');
% 
% 	Location = strcat('/main/DE/keypointdata/',modelData{k},'/SignalNoise/ShortRun/Std_');
% 
%     
% 	FileLocation = strcat(ProjectRoot,Location,'0.1','/');
%      
%     FileName = strcat('NMSKeypoint','.mat');
%     load(fullfile(FileLocation, FileName),'NMSKeypoint')
% 
%   
%     LevelTable = tabulate(NMSKeypoint.Level);
%     
%     KeypointLevelCounts(LevelTable(:,1),k) = LevelTable(:,2);
%     
%     KeypointLevels(1:3000,k) = (ScaleParameter) ./ PointCloud.Resolution;
% 
%     plot(KeypointLevels(:,k), KeypointLevelCounts(:,k)./sum(KeypointLevelCounts(:,k)), '-', 'linewidth', 7)
%     
% end
% legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'Orientation','Horizontal','FontSize',45) 
% xlim([0, 5])
% xlabel({'$\sigma_n / \bar{e}$'},'interpreter','Latex')
% ylabel('% of Keypoints')
% ax = gca;
% ax.XAxis.FontSize = 50;
% ax.YAxis.FontSize = 50;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
%     
%     
% 
% KeypointLevelCounts = zeros(3000,modelDataLength);
% % KeypointScaleCounts = cell(modelDataLength,1);
% markers = {'+-','o-','d-','s-','x-'};
% a = [10,5,1];
% figure(5)
% hold on
% 
% for k = 1 : modelDataLength
%     
%     
%     Model = strcat(modelData{k},'/', ModelNames{k})
%     FileLocationModel = strcat(ProjectRoot,'/models/object/');
%     FileNameModelPly = strcat(Model,'.ply');
%     
%     [PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
%     PointCloud = findMeshResolution(PointCloud, 'Model');
%     PointCloud
%     
%     tau = PointCloud.Resolution / 10;   
%     
%     ScaleParameter = findScaleParamter(tau, 1, 3000, 'Laplacian', 'Natural');
% 
% 	FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',modelData{k},'/SignalNoise/ShortRun/');
% 
%     FileName = strcat('NMSKeypoint','.mat');
%     load(fullfile(FileLocation, FileName),'NMSKeypoint')
% 
%     
%     LevelTable = tabulate(NMSKeypoint.Level);
%     
%     KeypointLevelCounts(LevelTable(:,1),k) = LevelTable(:,2);
%     
%     KeypointLevels(1:3000,k) = (ScaleParameter) ./ PointCloud.Resolution;
% 
%     plot(KeypointLevels(:,k), KeypointLevelCounts(:,k), '.-')
%     
% end
% legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'Orientation','Horizontal') 
% xlim([0, 5])
% xlabel({'$n \tau / \bar{e}$'},'interpreter','Latex')
% xlabel({'$\sigma_i / \bar{e}$'},'interpreter','Latex')
% ylabel('Count')
% ax = gca;
% ax.XAxis.FontSize = 25;
% ax.YAxis.FontSize = 25;
% 
% 
% 
% for kk = 1 : 3 %modelDataLength
%     plot(KeypointLevels(:,kk)*PointCloud.Resolution / (PointCloud.Resolution/a(kk)), KeypointLevelCounts(:,kk), '.-')%, 'Linewidth', 7, 'MarkerSize', 20)
% %     plot(1:3000, KeypointLevelCounts(:,kk), '.-')%, 'Linewidth', 7, 'MarkerSize', 20)
% end
% % legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'Orientation','Horizontal') 
% 
% xlim([0,150])
% xlabel({'$\frac{Level}{\bar{e}}$'},'interpreter','Latex')
% ylabel('Count')
% ax = gca;
% ax.XAxis.FontSize = 25;
% ax.YAxis.FontSize = 25;
% % legend({'Bunny','Armadillo','Buddha','Dragon','Itokawa'},'FontSize',45,'Orientation','Horizontal')        
% 
% 
% 
% 
% % Find percentage of model size for a given number of steps
% n = 100;
% tauFraction = 1/10;
% for k = 1 : modelDataLength
%     
%     Model = strcat(modelData{k},'/', ModelNames{k})
%     FileLocationModel = strcat(ProjectRoot,'/models/object/');
%     FileNameModelPly = strcat(Model,'.ply');
%     
%     [PointCloud.Location, PointCloud.Face] = read_ply( fullfile( FileLocationModel, FileNameModelPly ) );
%     PointCloud = findMeshResolution(PointCloud, 'Model');
%     PointCloud.Resolution
%     
% %     100*sqrt(2*n*PointCloud.Resolution*tauFraction)./range(PointCloud.Location)
% end







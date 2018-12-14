% Andrew Rhodes
% ASEL
% December 2018


close all
clear
clc

global ProjectRoot;

UseNMS = 0;
UseLocalReso = 1;
% t_scale = 1/sqrt(2);
% level_min = 1;
% t_range = 1/2;
alpha = 1;
NumSteps =2000;

TitleNames = {'Mesh Geo.', 'Mesh Euc.', 'Umbrella', 'Cotangent'};

% LBO = {'Mesh', 'Mesh_euc', 'Umbrella', 'Cotangent'};
LBO = {'Mesh_psp_opt_rho6', 'Mesh_psp_opt_rho6_euc', 'Umbrella_psp_opt', 'Cotangent_psp_opt'};
NoiseType = 'Signal';
RunType = '_te8';
modelFolder = {'bunny','armadillo','buddha','dragon','itokawa'}; % 
modelNames = {'Bunny_e1','Armadillo_e1_100000','Buddha_e1_50000', 'Dragon_e1_50000', 'Itokawa_e1_80000'};
 
 
LBOLength = length(LBO);
modelFolderLength = length(modelFolder);
modelNamesLength = length(modelNames);

RhodesWDLocation = '/media/andrew/WDRhodes/diffusiondata/';

% setTau = @(e_bar) e_bar/8;
setTau = @(e_bar) e_bar^(2/5);

LineStyles = {'b-', 'g-', 'r-', 'k-', 'c-'};


for kk = 1 : LBOLength
    
    
    
    figure
    title( TitleNames{kk})
    hold on
%     xlabel({'$\sigma / \bar{e}$'}, 'interpreter','Latex')
    xlabel({'$\sigma$'}, 'interpreter','Latex')
    ylabel('Percentage of Max Bin')
    ax = gca;
    ax.XAxis.FontSize = 50;
    ax.XLabel.FontSize = 75;
    ax.YAxis.FontSize = 50;
    ax.Title.FontSize = 50;
    
    
    for k = 1 : modelFolderLength
        
        
        
        % load PointCloud
        ModelToLoad = strcat(ProjectRoot,'/models/object/',modelFolder{k},'/',modelNames{k},'.ply');
        [PointCloud.Location, PointCloud.Face, PointCloud.Normal, PointCloud.Signal]...
            = read_ply_all_elements( ModelToLoad );
        PointCloud.LocationCount = size(PointCloud.Location,1);
        PointCloud.FaceCount = size(PointCloud.Face,1);
        PointCloud = findMeshResolution(PointCloud, 'Model');
%         PointCloud
        
        % load neighbors
        FileLocationNeighbors = strcat(RhodesWDLocation,modelFolder{k},'/neighbors/');
        FileNameNeighbors = strcat(modelNames{k},'_Neighbors.mat');
        load(strcat( FileLocationNeighbors, FileNameNeighbors), 'Neighbors');
        PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);
        
        %%%%%%%
        tau = setTau(PointCloud.Resolution);
        %%%%%%%
        
%         KeypointName = strcat(RhodesWDLocation,modelFolder{k},'/',NoiseType,'Noise/',LBO{kk},RunType,'/');
        KeypointName = strcat(RhodesWDLocation,modelFolder{k},'/',NoiseType,'Noise/',LBO{kk},'/');
        if UseNMS
            load( strcat(KeypointName, 'NMSKeypoint.mat') )
            Keypoint = NMSKeypoint;
            clear NMSKeypoint
        else
            load( strcat(KeypointName, 'Keypoint.mat') )
        end
        
        if UseLocalReso
            ScalePerResolution = Keypoint.Scale ./ PointCloud.ResolutionLocal(Keypoint.LocationIndex);
        else
            ScalePerResolution = Keypoint.Scale ./ PointCloud.Resolution;
        end
        
        ScaleParameter = findScaleParameter(tau, alpha, NumSteps, 'Laplacian', 'Natural');
        
        
%         [HistHeight, HistEdge, HistBin] = histcounts(Keypoint.Level, NumSteps, 'BinWidth', 1);
        
        LevelTable = tabulate(Keypoint.Level);
        
% % % %         figure
%         plot(ScaleParameter(LevelTable(:,1))./PointCloud.Resolution, ...
%             LevelTable(:,2)./max(LevelTable(:,2)), LineStyles{k}, 'Linewidth', 4)
        plot(ScaleParameter(LevelTable(:,1)), ...
            LevelTable(:,2)./max(LevelTable(:,2)), LineStyles{k}, 'Linewidth', 4)
        drawnow
        hold on
%         pause
        xlim([0 10])
        
        
        %
        %
        %     [HistHeight, HistEdge, HistBin] = histcounts(ScalePerResolution, length(ScalePerResolution));
        %
        %
        %     histogram(ScalePerResolution, length(ScalePerResolution));
        %
        %     figure
        %     plot(1:length(HistHeight), HistHeight)
        %     ax = gca;
        %     ax.XTickLabel = (1:length(ScalePerResolution)) * PointCloud.Resolution / 8;
        %
        %
        
        clear Keypoint PointCloud LevelTable tau
    
    end
    
    
    
end






















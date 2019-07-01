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
% NumSteps = 2000;

TitleNames = {'Mesh Geo.', 'Mesh Euc.', 'Umbrella', 'Cotangent'};

% LBO = {'Mesh_te8', 'Mesh_euc_te8', 'Umbrella_te8', 'Cotangent_te8'};
LBO = {'Mesh_rho4_ddr_geo', 'Mesh_rho4_ddr_euc', 'Umb2', 'Cot'};
% LBO = {'Mesh_rho4_ddr_geo_sph','Mesh_rho4_ddr_euc_sph' 'Umb_sph', 'Cot_sph'};
% LBO = {'Mesh_psp_opt_rho6', 'Mesh_psp_opt_rho6_euc', 'Umbrella_psp_opt', 'Cotangent_psp_opt'};
NoiseType = 'Vertex'; % 'Signal', 'Vertex'
modelFolder = {'bunny','armadillo','buddha','dragon','itokawa'}; % 
modelNames = {'Bunny_e1','Armadillo_e1_100000','Buddha_e1_50000', 'Dragon_e1_50000', 'Itokawa_e1_80000'};
 
 
LBOLength = length(LBO);
modelFolderLength = length(modelFolder);
modelNamesLength = length(modelNames);

RhodesWDLocation = '/media/andrew/WDRhodes/diffusiondata/';

% setTau = @(e_bar) e_bar/8;
% setTau = @(e_bar) 0.25;
konstant = 1.2;
l_ebar = 80;
set_ltau = @(e_bar) 4*e_bar;
set_t0 = @(e_bar) e_bar/4 ;
set_NumSteps = @(e_bar, t_0) ceil(log((l_ebar*e_bar)^2 / (2*alpha*t_0)) / (2*log(konstant)));

LineStyles = {'b-', 'g-', 'r-', 'k-', 'c-'};


for kk = 1 : LBOLength
    
    
    
    figure
    title( TitleNames{kk})
    hold on
    if kk == 3
        xlabel({'$\sigma$'}, 'interpreter','Latex')
    elseif kk == 4
        xlabel({'$\sigma / \sqrt{\bar{e}}$'}, 'interpreter','Latex')
    else
        xlabel({'$\sigma / \bar{e}$'}, 'interpreter','Latex')
    end
    ylabel('Percentage of Max Bin')
    ax = gca;
    ax.XAxis.FontSize = 50;
    ax.XLabel.FontSize = 75;
    ax.YAxis.FontSize = 50;
    ax.Title.FontSize = 50;
    
    x = [];
    for k = 1 : modelFolderLength
        
        
        
        % load PointCloud
        ModelToLoad = strcat(ProjectRoot,'/models/object/',modelFolder{k},'/',modelNames{k},'.ply');
        [PointCloud.Location, PointCloud.Face, PointCloud.Normal, PointCloud.Signal]...
            = read_ply_all_elements( ModelToLoad );
        PointCloud.LocationCount = size(PointCloud.Location,1);
        PointCloud.FaceCount = size(PointCloud.Face,1);
        PointCloud = findMeshResolution(PointCloud, 'Model');
        
        l_tau = set_ltau(PointCloud.Resolution);
        if kk == 3 % umbrella
            t_0 = set_t0(1);
            NumSteps = set_NumSteps(1, set_t0(1));
        else                    
            t_0 = set_t0(PointCloud.Resolution);
            NumSteps = set_NumSteps(PointCloud.Resolution, set_t0(PointCloud.Resolution));
        end
       

%         PointCloud
        
        % load neighbors
        FileLocationNeighbors = strcat(RhodesWDLocation,modelFolder{k},'/neighbors/');
        FileNameNeighbors = strcat(modelNames{k},'_Neighbors.mat');
        load(strcat( FileLocationNeighbors, FileNameNeighbors), 'Neighbors');
        PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);
        

        
        KeypointName = strcat(ProjectRoot,'/main/DE/keypointdata/',modelFolder{k},'/',NoiseType,'Noise/',LBO{kk},'/');
%         KeypointName = strcat(RhodesWDLocation,modelFolder{k},'/',NoiseType,'Noise/',LBO{kk},'/');
        if UseNMS
            load( strcat(KeypointName, 'NMSKeypoint.mat') )
            Keypoint = NMSKeypoint;
            clear NMSKeypoint
        else
            load( strcat(KeypointName, 'Keypoint.mat') )
        end
        
%         ZeroLogic = Keypoint.Scale < 3*PointCloud.Resolution;
% 
%         FNames = fieldnames(Keypoint);
%         for jj = 1 : length(FNames)
%             if strcmpi(FNames{jj},'Count')
%                 Keypoint = rmfield(Keypoint, FNames{jj});
%             else
%                 Keypoint.(FNames{jj})(ZeroLogic,:) = [];
%             end
%         end
%         Keypoint.Count = length(Keypoint.LocationIndex)


       
%         if UseLocalReso
%             ScalePerResolution = Keypoint.Scale ./ PointCloud.ResolutionLocal(Keypoint.LocationIndex);
%         else
%             ScalePerResolution = Keypoint.Scale ./ PointCloud.Resolution;
%         end
    


        [tn, tau_n, sigma_n] = findScaleStep(konstant, t_0, alpha, NumSteps);
%         ScaleParameter = findScaleParameter(tau, alpha, NumSteps, 'Laplacian', 'Natural');  
%         ScaleParameter = sqrt((0:NumSteps-1)'  * alpha * tau);
        
%         [HistHeight, HistEdge, HistBin] = histcounts(Keypoint.Level, NumSteps, 'BinWidth', 1);
        
        LevelTable = tabulate(Keypoint.Level);
        
        [a,b] = max(LevelTable(:,2));

        if kk == 3 % Umbrella LBO
            x = [x; sigma_n(LevelTable(b,1))];
            plot(sigma_n(LevelTable(:,1)), ...
                LevelTable(:,2)./max(LevelTable(:,2)), LineStyles{k}, 'Linewidth', 4)
        elseif kk == 4 % Cotangent LBO
            x = [x; sigma_n(LevelTable(b,1))./sqrt(PointCloud.Resolution)];
            plot(sigma_n(LevelTable(:,1))./sqrt(PointCloud.Resolution), ...
                LevelTable(:,2)./max(LevelTable(:,2)), LineStyles{k}, 'Linewidth', 4)
        else % Mesh LBO
            x = [x; sigma_n(LevelTable(b,1))./PointCloud.Resolution];
            plot(sigma_n(LevelTable(:,1))./PointCloud.Resolution, ...
                LevelTable(:,2)./max(LevelTable(:,2)), LineStyles{k}, 'Linewidth', 4)
            
%             x = [x; ScaleParameter(LevelTable(b,1))];
%             plot(ScaleParameter(LevelTable(:,1))./PointCloud.Resolution, ...
%                 LevelTable(:,2)./max(LevelTable(:,2)), LineStyles{k}, 'Linewidth', 4)
        end
        
        drawnow
        hold on
%         xlim([0 5])
        
%         pause
        
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
    title([])
    xlim([min(sigma_n) 5])
    
    x'
    sprintf('%s: mean: %0.4f   std: %0.4f', TitleNames{kk}, mean(x), std(x))
     
%     pause
end






















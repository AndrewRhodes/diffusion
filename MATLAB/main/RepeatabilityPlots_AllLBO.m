

close all
clear
clc

global ProjectRoot; % Additional Paths

UseNMS = 1;
t_scale = 1/sqrt(2);
level_min = 3;
t_range = 1/2;
NoiseType = 'Signal'; % 'Signal', 'Vertex' %
LBO = {'Mesh_','Mesh_euc','Umbrella','Cotangent'};
% LBO = {'Mesh_psp_opt_rho6_geo_LBM1','Mesh_psp_opt_rho6_euc','Umbrella_psp_opt','Cotangent_psp_opt'};
RunType = '_te8'; % 'Umbrella', 'Cotangent', 'Mesh', 'Mesh_euc_te8', 'Mesh_te8' %
lines = 1;
ByModel = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Keypoint Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelData = {'bunny','armadillo','buddha','dragon','itokawa'}; %{'dragon'}; % {'bunny','armadillo','buddha','dragon','itokawa'}; % {'buddha','dragon'}; %
ModelNames = {'Bunny_e1','Armadillo_e1_100000','Buddha_e1_50000', 'Dragon_e1_50000', 'Itokawa_e1_80000'}; %{'Dragon_e1_50000'}; %   {'Bunny_e1','Armadillo_e1_100000','Buddha_e1_50000', 'Dragon_e1_50000', 'Itokawa_e1_80000'};  % {'Buddha_e1_50000','Dragon_e1_50000'}; %

modelDataLength = length(modelData);
LBOLength = length(LBO);
Error = cell(modelDataLength, modelDataLength, LBOLength);
NoMatch = cell(modelDataLength, modelDataLength, LBOLength);
MultipleMatch = cell(modelDataLength, modelDataLength);
NoiseVec = [0.1, 0.2, 0.3, 0.4, 0.5];
NoiseVecLength = length(NoiseVec);





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
    
    load(strcat(FileLocationNeighbors, ModelNames{k},'_Neighbors.mat'), 'Neighbors')
    PointCloud = findLocalResolution(PointCloud, Neighbors.Connect);
    PointCloud = findMeshNormals(PointCloud);
    
    
    for kkk = 1 : LBOLength
        
        FileLocation = strcat(ProjectRoot,'/main/DE/keypointdata/',modelData{k},'/',NoiseType,'Noise/',LBO{kkk},RunType,'/');
%         FileLocation = strcat('/media/andrew/WDRhodes/diffusiondata/',modelData{k},'/',NoiseType,'Noise/',LBO{kkk},RunType,'/');
%         FileLocation = strcat('/media/andrew/WDRhodes/diffusiondata/',modelData{k},'/',NoiseType,'Noise/',LBO{kkk},'/');
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

            if strcmpi(NoiseType,'Vertex')
                
                [Error{k, kk, kkk}, NoMatch{k, kk, kkk}] = findKeypointErrorVertexNoise(PointCloud, NoiseVec(kk), ...
                    FileLocation, NoiseFileLocation, NumIter, UseNMS, t_scale, level_min, t_range);
                
            else

                [Error{k, kk, kkk}, NoMatch{k, kk, kkk}] = findKeypointErrorNew(PointCloud, FileLocation,...
                    NoiseFileLocation, NumIter, UseNMS, t_scale, level_min, t_range);

            end

        end
    end
    
    
    clear PointCloud
    
    
end



positions = [1, 1.2, 1.4, 1.6, 1.8];
% LineStyles = {'k-', 'r--', 'b:'};
LineStyles = {'k-', 'g--', 'r:', 'b-.'};
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
for jjj = 1 : LBOLength
    for jj = 1 : modelDataLength
        cp = positions + (jj - 1)*2;
        Data = zeros(NoiseVecLength,1);
        for j = 1 : NoiseVecLength
            Data(j) = median(Error{jj, j, jjj}.ScaleRepeat);
        end
        p(jjj) = plot(cp, Data, LineStyles{jjj}, 'Linewidth', 4);
    end
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
% legend(p, {'Mesh Geo.','Umb.','Cot.'}, 'FontSize', 30)
legend(p, {'Mesh Geo.','Mesh Euc.','Umb.','Cot.'}, 'FontSize', 50, 'Orientation','Horizontal', 'Location', 'northoutside')
% legend(p, {'Mesh Geo.','Mesh Euc.','Umb.','Cot.'}, 'FontSize', 50, 'Orientation','Horizontal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
for jjj = 1 : LBOLength
    for jj = 1 : modelDataLength
        cp = positions + (jj - 1)*2;
        Data = zeros(NoiseVecLength,1);
        for j = 1 : NoiseVecLength
            Data(j) = median(Error{jj, j, jjj}.RelativeRepeat);
        end
        p(jjj) = plot(cp, Data, LineStyles{jjj}, 'Linewidth', 4);
    end
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
% legend(p, {'Mesh Geo.','Umb.','Cot.'}, 'FontSize', 30)
legend(p, {'Mesh Geo.','Mesh Euc.','Umb.','Cot.'}, 'FontSize', 50, 'Orientation','Horizontal', 'Location', 'northoutside')
% legend(p, {'Mesh Geo.','Mesh Euc.','Umb.','Cot.'}, 'FontSize', 50, 'Orientation','Horizontal')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
hold on
for jjj = 1 : LBOLength
    for jj = 1 : modelDataLength
        cp = positions + (jj - 1)*2;
        Data = zeros(NoiseVecLength,1);
        for j = 1 : NoiseVecLength
            Data(j) = median(Error{jj, j, jjj}.Count);
        end
        p(jjj) = plot(cp, Data, LineStyles{jjj}, 'Linewidth', 4);
        ylims(jj,:) = get(gca, 'YLim');
    end
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
% legend(p, {'Mesh Geo.','Umb.','Cot.'}, 'FontSize', 30)
legend(p, {'Mesh Geo.','Mesh Euc.','Umb.','Cot.'}, 'FontSize', 50, 'Orientation','Horizontal', 'Location', 'northoutside')
% legend(p, {'Mesh Geo.','Mesh Euc.','Umb.','Cot.'}, 'FontSize', 50, 'Orientation','Horizontal')






